# rmetis 技术规范文档

> 版本: 0.1.1 | 状态: 草稿 | 日期: 2026-04-03

## 1. 项目概述

rmetis 是一个纯 Rust 实现的图划分库，功能对标 METIS 5.1.x 和 ParMETIS 4.x。目标：

- 完全纯 Rust，无 C/C++ 依赖
- 可编译为 WebAssembly (`wasm32-unknown-unknown`)
- 兼容 METIS 的核心 API 语义
- 支持 `no_std`（可选特性）
- MPI 级别分布式并行（ParMETIS 风格并行 ncuts、分布式粗化与细化）

---

## 2. 图划分算法规范

### 2.1 多层图划分框架（Multilevel Framework）

所有划分算法均基于三阶段多层框架：

```
原始图 G₀
    ↓  [粗化阶段] 匹配 + 合并顶点
   G₁ → G₂ → ... → Gₘ  (|Gₘ| ≤ 粗化阈值)
    ↓  [初始划分] 对 Gₘ 直接划分
   Partition(Gₘ)
    ↓  [细化阶段] 逐层投影 + 局部优化
   Partition(Gₘ₋₁) → ... → Partition(G₀)
```

**粗化停止条件：**
- `|V(Gᵢ)| ≤ max(20 × k, coarsen_limit)`，其中 `k` 为目标分区数
- 连续两层粗化率 < 5%（防止无效粗化循环）

### 2.2 粗化算法（Coarsening）

#### 2.2.1 随机匹配（Random Matching, RM）

```
算法 RM-Coarsen(G):
  visited = [false; nvtxs]
  matching = [-1; nvtxs]
  permutation = random_shuffle(0..nvtxs)
  
  for v in permutation:
    if visited[v]: continue
    visited[v] = true
    for u in neighbors(v):
      if not visited[u]:
        matching[v] = u
        matching[u] = v
        visited[u] = true
        break
  
  return build_coarse_graph(G, matching)
```

#### 2.2.2 排序重边匹配（Sorted Heavy Edge Matching, SHEM）

```
算法 SHEM-Coarsen(G):
  visited = [false; nvtxs]
  matching = [-1; nvtxs]
  
  // 预先采样每个顶点的随机平局打破值（在排序前一次性生成）
  // 注意：不能在 sort_by_key 闭包内调用 rng.gen()，
  // 因为 sort_by_key 对同一元素会多次调用闭包，违反全序关系。
  tiebreakers = [rng.gen::<u32>(); nvtxs]
  permutation = sort_by_key(0..nvtxs, key = (degree(v), tiebreaker[v]))
  
  for v in permutation:
    if visited[v]: continue
    visited[v] = true
    best_u = -1
    best_weight = -∞
    for (u, w) in weighted_neighbors(v):
      if not visited[u] and w > best_weight:
        best_u = u; best_weight = w
    if best_u != -1:
      matching[v] = best_u
      matching[best_u] = v
      visited[best_u] = true
  
  return build_coarse_graph(G, matching)
```

#### 2.2.3 粗化图构造

```
算法 build_coarse_graph(G, matching):
  // 为每个超顶点分配 ID
  cmap = [-1; nvtxs]
  cnvtxs = 0
  for v in 0..nvtxs:
    if cmap[v] == -1:
      cmap[v] = cnvtxs
      if matching[v] != -1:
        cmap[matching[v]] = cnvtxs
      cnvtxs++
  
  // 合并顶点权重
  cvwgt[u] = vwgt[v1] + vwgt[v2]  // v1, v2 匹配到超顶点 u
  
  // 合并边（去重+权重累加）
  for each super-vertex u:
    for each original v in u:
      for (w, ew) in neighbors(v):
        if cmap[w] != u:  // 跨超顶点的边
          cadjncy[u] += (cmap[w], ew)  // 累加重边权重
```

### 2.3 初始划分（Initial Partitioning）

#### 2.3.1 图增长法（Graph Growing, GG）

```
算法 GG-Partition(G, k):
  part = [-1; nvtxs]
  for p in 0..k:
    seed = random_unassigned_vertex()
    queue = BFS_queue([seed])
    target_weight = total_weight / k
    while queue not empty and part_weight[p] < target_weight:
      v = queue.pop()
      if part[v] == -1:
        part[v] = p
        part_weight[p] += vwgt[v]
        queue.extend(unassigned_neighbors(v))
  
  // 处理剩余未分配顶点
  assign_remainder_greedy(part)
  return part
```

#### 2.3.2 随机初始划分（Random）

均匀随机分配顶点到各分区，然后执行细化。

#### 2.3.3 递归二分初始划分（Metis-RB）

仅用于 k-way 算法的初始阶段，递归二分直到达到目标分区数。

### 2.4 细化算法（Refinement）

#### 2.4.1 FM（Fiduccia-Mattheyses）算法

**边割最小化版本（2-way）：**

```
算法 FM-Refine-2way(G, part, niter):
  for iter in 0..niter:
    gain = compute_initial_gains()    // gain[v] = 切割减少量
    pq = MaxPriorityQueue(boundary_vertices, gain)
    locked = [false; nvtxs]
    best_cut = current_cut
    best_state = clone(part)
    cumulative_gain = 0
    moves = []
    
    while pq not empty:
      v = pq.pop_max()
      if locked[v]: continue
      if move_improves_balance(v):
        do_move(v, part)
        locked[v] = true
        cumulative_gain += gain[v]
        moves.push(v)
        update_gains(neighbors(v), pq)  // 更新邻居的增益
        if current_cut < best_cut:
          best_cut = current_cut
          best_state = clone(part)
    
    part = best_state  // 回滚到最佳状态
    if no_improvement: break
```

**增益计算：**
```
gain[v] = (edges from v to other_part) - (edges from v to same_part)
        = cut_edges(v) - internal_edges(v)
```

#### 2.4.2 贪心细化（Greedy Boundary Refinement）

贪心细化分两个阶段：

**阶段一：平衡修复 pass（balance-repair）**

在切割优化循环之前先运行一次，将顶点从超重分区移出，即使移动不改善切割也允许执行：

```
算法 Balance-Repair(G, part, nparts, tpwgts, ubvec):
  for v in boundary_vertices:
    if source_partition(v) 不超重: continue
    best_p = argmax_gain(v, dest in feasible_partitions)
    // 可行性：移入后目标分区不超过 tpwgts[p] × ubvec 阈值
    if best_p exists:
      move(v, best_p)
      update partition weights
```

**阶段二：切割优化 passes（niter 次迭代）**

```
算法 Greedy-Refine(G, part, niter):
  for iter in 0..niter:
    improved = false
    for v in boundary_vertices:
      best_part = argmax_gain(v, dest satisfying balance constraint)
      if gain(v, best_part) > 0:
        move(v, best_part)
        improved = true
    if not improved: break
```

#### 2.4.3 K-way FM 细化

扩展 FM 算法至 k-way 场景：

```
gain[v][p] = (edges from v to partition p) - (edges from v to current_part(v))
```

每次移动选择 `argmax_p gain[v][p]`，约束余量平衡。

### 2.5 递归二分（Recursive Bisection）

```
算法 RB-Partition(G, k, tpwgts):
  if k == 1:
    return all_same_partition(G)
  
  k1 = k / 2; k2 = k - k1
  target = [sum(tpwgts[0..k1]), sum(tpwgts[k1..k])]
  
  part_2 = multilevel_bisect(G, target)
  
  G1, G2 = split_graph(G, part_2)
  
  part1 = RB-Partition(G1, k1, tpwgts[0..k1])
  part2 = RB-Partition(G2, k2, tpwgts[k1..k])
  
  return merge_partitions(part1, part2, k1)
```

### 2.6 嵌套剖分（Nested Dissection for NodeND）

```
算法 ND-Order(G, perm, iperm):
  if |G| ≤ nd_threshold:
    // 直接按任意顺序排列小图
    assign_sequential(G, perm, iperm)
    return
  
  separator = find_vertex_separator(G)  // 找分隔符
  G_left, G_right = G - separator
  
  // 递归处理子图（先排子图，再排分隔符）
  ND-Order(G_left, perm, iperm)
  ND-Order(G_right, perm, iperm)
  assign_next(separator, perm, iperm)   // 分隔符排在最后
```

### 2.7 MPI 并行算法

rmetis 通过 `Comm` trait 抽象 MPI 通信，支持三层并行：并行 ncuts、分布式粗化、分布式贪心细化。

#### 2.7.1 并行 ncuts（Parallel ncuts）

当 `comm.size() > 1` 时，将 `ncuts` 次独立尝试分发到各 MPI 进程：

```
算法 Parallel-ncuts(graph, nparts, ncuts, nprocs):
  // 进程 r 处理尝试编号 r, r+P, r+2P, ...
  for attempt in (rank .. ncuts).step_by(nprocs):
    seed = base_seed + attempt
    local_result = single_kway_attempt(graph, nparts, seed)
  
  // 全局归约：找最小边割
  global_best_cut = all_reduce_min(local_best_cut)
  
  // 持有最优结果的进程将其 gather 到 rank 0，再 broadcast 给所有进程
  winning_part = gather_then_broadcast(local_part if local_cut == global_best_cut)
```

#### 2.7.2 分布式粗化（Distributed Coarsening）

每个进程只为其"拥有"的顶点（`v % nprocs == rank`）计算匹配，然后通过 `all_gather` 交换：

```
算法 Distributed-Coarsen(G, nprocs, rank):
  // 独立计算本进程负责顶点的匹配建议
  local_proposals = []
  for v in 0..nvtxs where v % nprocs == rank:
    best_u = argmax_weight_neighbor(v, unvisited)
    if best_u found:
      local_proposals.push((v, best_u))
  
  // 交换所有进程的建议
  all_proposals = all_gather(local_proposals)
  
  // 冲突消解：按进程顺序处理，最小 rank 的建议优先
  matching = identity_matching(nvtxs)
  claimed = [false; nvtxs]
  for rank_proposals in all_proposals:  // 按 rank 顺序
    for (v, u) in rank_proposals:
      if not claimed[v] and not claimed[u]:
        matching[v] = u; matching[u] = v
        claimed[v] = true; claimed[u] = true
  
  // 所有进程得到相同的 matching，本地构造粗化图
  return build_coarse_graph(G, matching)
```

#### 2.7.3 分布式贪心细化（Distributed Greedy Refinement）

边界顶点按 `v % nprocs == rank` 分配给各进程。每轮两个阶段（平衡修复 + 切割优化），均通过 `all_gather` 原子交换移动：

```
算法 Distributed-Greedy-Refine(G, part, nprocs, rank, niter):
  // 平衡修复 pass
  local_moves = propose_balance_repair_moves(rank)
  all_moves = all_gather(local_moves)
  apply_moves(all_moves)
  sync_partition_weights()
  
  // 切割优化 passes
  for iter in 0..niter:
    local_moves = propose_cut_improving_moves(rank)
    all_moves = all_gather(local_moves)
    any_move = apply_moves(all_moves)
    if not any_move: break
    sync_partition_weights()
```

---

## 3. 数据结构规范

### 3.1 核心图结构（CSR 格式）

```rust
/// 压缩稀疏行图表示
pub struct Graph {
    /// 顶点数
    pub nvtxs: usize,
    /// 边数（无向图中每条边存储两次，此处为实际边数）
    pub nedges: usize,
    /// 邻接索引数组，长度 nvtxs + 1
    /// xadj[i]..xadj[i+1] 为顶点 i 的邻居在 adjncy 中的范围
    pub xadj: Vec<Idx>,
    /// 邻接数组，长度 2 * nedges（无向图）
    pub adjncy: Vec<Idx>,
    /// 顶点权重（多约束：长度 nvtxs * ncon），None 表示均匀权重 1
    pub vwgt: Option<Vec<Idx>>,
    /// 边权重，长度 2 * nedges，None 表示均匀权重 1
    pub adjwgt: Option<Vec<Idx>>,
    /// 顶点尺寸（通信体积优化用），None 表示均匀尺寸 1
    pub vsize: Option<Vec<Idx>>,
    /// 约束数量（多约束划分时 > 1）
    pub ncon: usize,
}
```

### 3.2 多层图层次结构

```rust
pub struct CoarseGraph {
    pub graph: Graph,
    /// 原始图顶点到此粗化图超顶点的映射
    pub cmap: Vec<Idx>,
    /// 原始层（用于回溯投影）
    pub finer: Option<Box<CoarseGraph>>,
}
```

### 3.3 划分结果

```rust
pub struct PartitionResult {
    /// 每个顶点所属分区 ID (0..nparts-1)
    pub part: Vec<Idx>,
    /// 目标函数值（边割数或通信体积）
    pub objval: Idx,
}
```

### 3.4 嵌套剖分结果

```rust
pub struct NDResult {
    /// 填充减少置换：perm[i] = 原始顶点 i 在新排序中的位置
    pub perm: Vec<Idx>,
    /// 逆置换：iperm[i] = 新排序位置 i 对应的原始顶点
    pub iperm: Vec<Idx>,
}
```

### 3.5 类型定义

```rust
/// 索引类型，默认 i32 以兼容 METIS C API
pub type Idx = i32;
/// 实数类型，用于权重比例
pub type Real = f32;
```

### 3.6 MPI 通信抽象（`comm.rs`）

```rust
/// MPI 通信抽象 trait，统一屏蔽 SingleComm / jsmpi / mpi 后端的差异。
/// 所有方法在单进程（SingleComm）下均为恒等操作或本地操作，无网络开销。
pub trait Comm: Send + Sync {
    /// 当前进程编号（0-based）
    fn rank(&self) -> i32;
    /// 总进程数
    fn size(&self) -> i32;
    /// 全局归约：取所有进程 local 中的最小值
    fn all_reduce_min_i32(&self, local: Idx) -> Idx;
    /// 全局归约（原位）：各进程对应元素求和，结果写回 out
    fn all_reduce_sum_i32_slice(&self, local: &[Idx], out: &mut [Idx]);
    /// 广播：rank root 将 data 广播给所有进程
    fn broadcast_i32_vec(&self, root: i32, data: &mut Vec<Idx>);
    /// 收集：各进程将 local 发往 root，root 返回 Vec<Vec<Idx>>
    fn gather_i32_vec(&self, root: i32, local: &[Idx]) -> Vec<Vec<Idx>>;
    /// 全量收集：每个进程收到所有进程的数据
    fn all_gather_i32_vec(&self, local: &[Idx]) -> Vec<Vec<Idx>>;
    /// 同步屏障
    fn barrier(&self);
}

/// 无操作单进程实现（不依赖 MPI）
pub struct SingleComm;

/// 创建世界通信器（自动选择后端）：
/// - wasm32 目标：使用 jsmpi 后端
/// - 其他目标：使用 mpi crate 后端（OpenMPI / MPICH）
/// - 若仅一个进程：返回 SingleComm（避免 MPI 初始化开销）
pub fn world() -> Box<dyn Comm>;
```

**后端选择规则：**

| 构建目标 | MPI 后端 | 说明 |
|---------|---------|------|
| `wasm32-*` | jsmpi（`vendor/jsmpi` 子模块） | Web Worker 间通信 |
| 其他目标，`size > 1` | `mpi` crate 0.8（OpenMPI / MPICH） | 原生 MPI |
| 其他目标，`size == 1` | `SingleComm`（no-op） | 无 MPI 开销 |

**注意**：jsmpi 原生只支持 sum-reduce，`all_reduce_min` 在 jsmpi 后端通过 `gather + local-min + broadcast` 实现。

---

## 4. 公开 API 规范

### 4.1 Rust 原生 API

```rust
/// 多层递归二分图划分
pub fn part_graph_recursive(
    graph: &Graph,
    nparts: usize,
    tpwgts: Option<&[Real]>,   // 长度 nparts * ncon，None 表示均匀
    ubvec: Option<&[Real]>,    // 长度 ncon，None 表示默认容差
    options: &Options,
) -> Result<PartitionResult, MetisError>;

/// 多层 k-way 图划分（串行，内部 comm = None）
pub fn part_graph_kway(
    graph: &Graph,
    nparts: usize,
    tpwgts: Option<&[Real]>,
    ubvec: Option<&[Real]>,
    options: &Options,
) -> Result<PartitionResult, MetisError>;

/// 多层 k-way 图划分（MPI 并行版本）
/// comm = Some(&world) 时启用并行 ncuts + 分布式细化
/// comm = None 时退化为串行（等同于 part_graph_kway）
pub fn partition_kway(
    graph: &Graph,
    nparts: usize,
    tpwgts: Option<&[Real]>,
    ubvec: Option<&[Real]>,
    options: &Options,
    comm: Option<&dyn Comm>,
) -> Result<PartitionResult, MetisError>;

/// 嵌套剖分（稀疏矩阵填充减少排序）
pub fn node_nd(
    graph: &Graph,
    options: &Options,
) -> Result<NDResult, MetisError>;
```

### 4.2 METIS C ABI 兼容层（`c-api` feature）

```c
/* 与 METIS 5.1.x 相同的函数签名 */
int RMETIS_PartGraphRecursive(
    idx_t *nvtxs, idx_t *ncon,
    idx_t *xadj, idx_t *adjncy,
    idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
    idx_t *nparts, real_t *tpwgts, real_t *ubvec,
    idx_t *options,
    idx_t *objval, idx_t *part
);

int RMETIS_PartGraphKway(
    idx_t *nvtxs, idx_t *ncon,
    idx_t *xadj, idx_t *adjncy,
    idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
    idx_t *nparts, real_t *tpwgts, real_t *ubvec,
    idx_t *options,
    idx_t *objval, idx_t *part
);

int RMETIS_NodeND(
    idx_t *nvtxs, idx_t *xadj, idx_t *adjncy,
    idx_t *vwgt, idx_t *options,
    idx_t *perm, idx_t *iperm
);

void RMETIS_SetDefaultOptions(idx_t *options);
```

### 4.3 WASM JavaScript API（`wasm` feature）

```typescript
// TypeScript 绑定（wasm-bindgen 生成）
interface RMetisGraph {
  nvtxs: number;
  xadj: Int32Array;
  adjncy: Int32Array;
  vwgt?: Int32Array;
  adjwgt?: Int32Array;
  ncon?: number;
}

interface PartitionOptions {
  ctype?: 'rm' | 'shem';          // 粗化类型
  iptype?: 'grow' | 'random';     // 初始划分类型
  rtype?: 'fm' | 'greedy';        // 细化类型
  objtype?: 'cut' | 'vol';        // 目标函数（仅 kway）
  ncuts?: number;                  // 尝试次数（取最优）
  niter?: number;                  // 每层细化迭代次数
  seed?: number;                   // 随机种子
  ufactor?: number;                // 不平衡容差 (1-1000)
  minconn?: boolean;               // 最小化分区连通度
  contig?: boolean;                // 强制连通分区
}

export function partGraphRecursive(
  graph: RMetisGraph,
  nparts: number,
  options?: PartitionOptions
): { part: Int32Array; objval: number };

export function partGraphKway(
  graph: RMetisGraph,
  nparts: number,
  options?: PartitionOptions
): { part: Int32Array; objval: number };

export function nodeNd(
  graph: RMetisGraph,
  options?: PartitionOptions
): { perm: Int32Array; iperm: Int32Array };
```

---

## 5. 配置选项规范

```rust
/// METIS_NOPTIONS = 40
pub struct Options {
    /// 粗化方案
    pub ctype: CType,
    /// 初始划分方案
    pub iptype: IpType,
    /// 细化方案
    pub rtype: RType,
    /// 目标函数类型（仅 k-way 使用）
    pub objtype: ObjType,
    /// 禁止 2-hop 匹配（0=启用，1=禁止）
    pub no2hop: bool,
    /// 划分尝试次数（取最优），范围 [1, 100]
    pub ncuts: usize,
    /// 每层细化迭代次数，范围 [1, 100]
    pub niter: usize,
    /// 随机种子（0=基于时间）
    pub seed: u64,
    /// 最小化分区最大连接度（仅 k-way）
    pub minconn: bool,
    /// 强制分区连通（仅 k-way）
    pub contig: bool,
    /// 不平衡容差 × 1000，范围 [1, 1000]，默认 1
    pub ufactor: u32,
    /// 顶点编号方式（0=0-based，1=1-based）
    pub numbering: u8,
    /// 调试级别（0=关闭，位掩码）
    pub dbglvl: u32,
}

#[derive(Default)]
pub enum CType { #[default] Shem, Rm }

#[derive(Default)]
pub enum IpType { #[default] Grow, Random, Edge, Node, MetisRb }

#[derive(Default)]
pub enum RType { #[default] Fm, Greedy, Sep2sided, Sep1sided }

#[derive(Default)]
pub enum ObjType { #[default] Cut, Vol }
```

---

## 6. 错误处理规范

```rust
#[derive(Debug, thiserror::Error)]
pub enum MetisError {
    #[error("输入错误: {0}")]
    Input(String),

    #[error("内存不足")]
    OutOfMemory,

    #[error("图不连通，且要求连通分区")]
    Disconnected,

    #[error("无法满足平衡约束: 最小不平衡={0:.4}, 要求={1:.4}")]
    UnbalancedPartition(f32, f32),

    #[error("内部错误: {0}")]
    Internal(String),
}

/// METIS 返回码（C ABI 层使用）
pub const METIS_OK: i32 = 1;
pub const METIS_ERROR_INPUT: i32 = -2;
pub const METIS_ERROR_MEMORY: i32 = -3;
pub const METIS_ERROR: i32 = -4;
```

---

## 7. 多约束划分规范

当 `ncon > 1` 时，每个顶点有多个权重维度，平衡约束变为向量形式：

```
对于每个约束 c 和每个分区 p:
  weight(p, c) / total_weight(c) ≤ tpwgts[p * ncon + c] × ubvec[c]
```

**细化时的多约束增益计算：**

移动顶点 v 从分区 s 到分区 t 的可行性：
```
for each constraint c:
  new_weight_s[c] = weight(s, c) - vwgt[v * ncon + c]
  new_weight_t[c] = weight(t, c) + vwgt[v * ncon + c]
  
  if new_weight_t[c] / total[c] > tpwgts[t * ncon + c] × ubvec[c]:
    move infeasible
```

---

## 8. WASM 约束与限制

### 8.1 内存模型

- WASM 线性内存默认 16MB，最大可扩展至 4GB
- 大图（> 100M 顶点）在浏览器环境中不适用
- 建议在 WASM 中处理 ≤ 10M 顶点的图

### 8.2 `wasm32-unknown-unknown` 要求

- 禁止使用 `std::thread`（可选：`wasm32-wasi` 支持线程）
- 禁止使用文件 I/O（除 WASI 外）
- 随机数：使用 `getrandom` crate（自动桥接浏览器 `crypto.getRandomValues`）
- 无 `std::time`：seed 由调用方传入

### 8.3 `no_std` 支持

```toml
[features]
default = ["std"]
std = ["thiserror/std", "rand/std"]
# no_std 时使用 alloc crate，随机数由外部注入
```

### 8.4 WASM 二进制尺寸优化目标

| 构建配置 | 目标大小 |
|---------|---------|
| Debug WASM | ≤ 5MB |
| Release WASM (wasm-opt -O2) | ≤ 500KB |
| Release WASM with LTO (wasm-opt -O3) | ≤ 300KB |

---

## 9. 性能基准规范

### 9.1 参考数据集

| 数据集 | 顶点数 | 边数 | 来源 |
|-------|--------|------|------|
| 4elt | 15,606 | 45,878 | METIS 测试集 |
| mdual | 258,569 | 513,562 | METIS 测试集 |
| copter2 | 9,167 | 49,152 | METIS 测试集 |
| auto | 448,695 | 3,314,611 | METIS 测试集 |

### 9.2 性能目标（原生 x86_64）

| 操作 | 图规模 | 目标时间 |
|------|--------|---------|
| PartGraphKway (k=8) | 100K 顶点 | ≤ 200ms |
| PartGraphKway (k=32) | 1M 顶点 | ≤ 5s |
| PartGraphRecursive (k=8) | 100K 顶点 | ≤ 300ms |
| NodeND | 100K 顶点 | ≤ 400ms |

### 9.3 划分质量目标

相比 METIS 5.1.x 参考实现，边割数差距 ≤ 5%（同等参数下）。

---

## 10. 测试规范

### 10.1 正确性测试

- **CSR 合法性验证**：对称性、无自环、索引范围
- **划分合法性**：所有顶点已分配、分区 ID 在范围内
- **平衡性验证**：每个分区权重不超过 `(1 + ufactor/1000) × target`
- **连通性验证**（`contig=true` 时）：BFS 验证每分区连通
- **NodeND 验证**：`perm` 和 `iperm` 互为逆置换
- **MPI Comm 验证**：`SingleComm` 与显式 `SingleComm` 结果一致、`all_gather` 返回单行、广播为恒等操作

### 10.2 回归测试

固定 `seed` 下，对标准数据集的边割数结果应确定性一致。

### 10.3 并行测试（`tests/parallel.rs`）

覆盖以下场景（均在 `SingleComm` 下运行，无需真实 MPI 环境）：

| 测试 | 验证内容 |
|------|---------|
| `single_comm_rank_and_size` | `rank()==0, size()==1` |
| `single_comm_all_reduce_min_is_identity` | 单进程 min = 原值 |
| `single_comm_all_reduce_sum_slice_is_copy` | 单进程 sum = 原切片 |
| `single_comm_gather_returns_single_row` | gather 返回 1 行 |
| `single_comm_all_gather_returns_single_row` | all_gather 返回 1 行 |
| `single_comm_broadcast_is_noop` | broadcast 不改变 data |
| `greedy_refine_with_single_comm_matches_serial` | comm 版与串行版结果相同 |
| `greedy_refine_with_comm_does_not_increase_cut` | 细化不增加切割数 |
| `coarsen_with_single_comm_produces_valid_hierarchy` | 粗化层次结构合法 |
| `partition_kway_with_explicit_single_comm_is_valid` | 显式 comm 划分合法 |
| `partition_kway_with_none_comm_matches_explicit_single_comm` | None 等价于 SingleComm |
| `partition_kway_parallel_ncuts_covers_all_attempts` | ncuts 尝试全部被执行 |

### 10.4 Fuzz 测试

使用 `cargo-fuzz` 对 CSR 输入进行模糊测试，验证不 panic、不越界。

---

## 11. Crate 结构

```
rmetis/
├── Cargo.toml
├── src/
│   ├── lib.rs              # 公开 API 导出
│   ├── types.rs            # Idx, Real, Graph, Options, 错误类型
│   ├── comm.rs             # MPI 通信抽象（Comm trait + SingleComm / jsmpi / mpi 后端）
│   ├── graph/
│   │   ├── mod.rs          # Graph 构造与验证
│   │   └── validate.rs     # CSR 合法性检查
│   ├── coarsen/
│   │   ├── mod.rs          # 粗化框架（coarsen + coarsen_with_comm）
│   │   ├── rm.rs           # 随机匹配
│   │   └── shem.rs         # SHEM（确定性平局打破）
│   ├── initial/
│   │   ├── mod.rs          # 初始划分调度
│   │   ├── grow.rs         # 图增长法
│   │   └── random.rs       # 随机划分
│   ├── refine/
│   │   ├── mod.rs          # 细化框架
│   │   ├── fm.rs           # FM 算法
│   │   ├── fm_kway.rs      # K-way FM
│   │   └── greedy.rs       # 贪心细化（greedy_refine + greedy_refine_with_comm）
│   ├── partition/
│   │   ├── mod.rs          # 划分入口
│   │   ├── recursive.rs    # 递归二分
│   │   ├── kway.rs         # K-way 划分（partition_kway，支持 comm 参数）
│   │   └── nd.rs           # 嵌套剖分
│   ├── separator/
│   │   └── mod.rs          # 顶点分隔符算法（NodeND 用）
│   └── ffi/
│       ├── c_api.rs        # C ABI 兼容层（feature: c-api）
│       └── wasm_api.rs     # WASM/JS 绑定（feature: wasm）
├── vendor/
│   └── jsmpi/              # git submodule — WASM 目标的 MPI 兼容层
├── tests/
│   ├── correctness.rs      # 算法正确性与回归测试
│   ├── parallel.rs         # MPI 抽象与并行算法测试（基于 SingleComm）
│   └── fixtures/           # 测试图数据文件
└── benches/
    └── partition.rs
```

---

## 12. Cargo.toml 依赖规范

```toml
[package]
name = "rmetis"
version = "0.1.0"
edition = "2021"
rust-version = "1.75"

[dependencies]
# 随机数（支持 wasm32）
rand = { version = "0.8", default-features = false, features = ["small_rng"] }
getrandom = { version = "0.2", features = ["js"], optional = true }

# 优先队列（FM 算法）
priority-queue = "1.3"

# 错误处理
thiserror = { version = "1", default-features = false }

# WASM 绑定（可选）
wasm-bindgen = { version = "0.2", optional = true }
serde = { version = "1", features = ["derive"], optional = true }
serde-wasm-bindgen = { version = "0.6", optional = true }

# MPI 后端：wasm32 目标使用 jsmpi（git submodule at vendor/jsmpi）
[target.'cfg(target_arch = "wasm32")'.dependencies]
jsmpi = { path = "vendor/jsmpi" }

# MPI 后端：非 wasm32 目标使用 mpi crate（wraps OpenMPI / MPICH）
[target.'cfg(not(target_arch = "wasm32"))'.dependencies]
mpi = "0.8"

[dev-dependencies]
criterion = "0.5"

[features]
default = ["std"]
std = ["rand/std", "rand/std_rng"]
c-api = []
wasm = ["wasm-bindgen", "serde", "serde-wasm-bindgen", "getrandom"]
parallel = ["rayon"]
rayon = ["dep:rayon"]

[lib]
crate-type = ["cdylib", "rlib"]

[[bench]]
name = "partition"
harness = false
```

---

## 附录 A：METIS 选项数组索引

```rust
pub const METIS_OPTION_PTYPE:     usize = 0;
pub const METIS_OPTION_OBJTYPE:   usize = 1;
pub const METIS_OPTION_CTYPE:     usize = 2;
pub const METIS_OPTION_IPTYPE:    usize = 3;
pub const METIS_OPTION_RTYPE:     usize = 4;
pub const METIS_OPTION_DBGLVL:    usize = 5;
pub const METIS_OPTION_NITER:     usize = 6;
pub const METIS_OPTION_NCUTS:     usize = 7;
pub const METIS_OPTION_SEED:      usize = 8;
pub const METIS_OPTION_NO2HOP:    usize = 9;
pub const METIS_OPTION_MINCONN:   usize = 10;
pub const METIS_OPTION_CONTIG:    usize = 11;
pub const METIS_OPTION_COMPRESS:  usize = 12;
pub const METIS_OPTION_CCORDER:   usize = 13;
pub const METIS_OPTION_PFACTOR:   usize = 14;
pub const METIS_OPTION_NSEPS:     usize = 15;
pub const METIS_OPTION_UFACTOR:   usize = 16;
pub const METIS_OPTION_NUMBERING: usize = 17;
pub const METIS_NOPTIONS:         usize = 40;
```

## 附录 B：参考文献

1. Karypis, G., & Kumar, V. (1998). A fast and high quality multilevel scheme for partitioning irregular graphs. *SIAM Journal on Scientific Computing*, 20(1), 359-392.
2. Karypis, G., & Kumar, V. (1998). Multilevel k-way partitioning scheme for irregular graphs. *Journal of Parallel and Distributed Computing*, 48(1), 96-129.
3. Karypis, G., Schloegel, K., & Kumar, V. (2003). *ParMETIS: Parallel Graph Partitioning and Sparse Matrix Ordering Library*. Technical Report, University of Minnesota.
4. Karypis, G., & Kumar, V. (1998). A parallel algorithm for multilevel graph partitioning and sparse matrix ordering. *Journal of Parallel and Distributed Computing*, 48(1), 71-95.
5. METIS 5.1.0 Manual: https://karypis.github.io/glaros/files/sw/metis/manual.pdf
