# rmetis 设计开发指南

> 面向 AI Agent 的开发指导文档 | 版本: 0.1.0 | 日期: 2026-04-02

## 1. 项目目标与约束

### 1.1 核心目标

| 目标 | 描述 |
|------|------|
| 算法完备性 | 实现 METIS 5.1.x 的全部三个核心函数：PartGraphRecursive、PartGraphKway、NodeND |
| WASM 兼容 | 编译目标 `wasm32-unknown-unknown`，无平台依赖代码 |
| API 兼容性 | 提供可选的 C ABI 层，使现有 METIS 用户零成本迁移 |
| 正确性优先 | 初期以正确性和算法一致性为主，性能优化在 v0.2 进行 |

### 1.2 硬性约束

- **禁止使用 unsafe**（除 FFI C ABI 层的必要边界）
- **禁止 C/C++ 依赖**：不能 FFI 调用原版 METIS
- **WASM 兼容**：所有核心逻辑必须在 `wasm32-unknown-unknown` 下通过编译
- **确定性**：相同输入 + 相同 seed，输出必须完全一致
- **Rust 最低版本**：1.75（稳定版 `async fn in trait` + `impl Trait in return position`）

---

## 2. 开发阶段规划

### Phase 1：基础设施（优先完成）

**目标**：可运行的最小框架，能对小图做简单二分。

- [ ] **P1.1** 建立 Cargo workspace，配置所有 feature flags
- [ ] **P1.2** 实现 `Graph` 结构体与 CSR 验证
- [ ] **P1.3** 实现随机匹配粗化（RM）
- [ ] **P1.4** 实现图增长初始划分（GG）
- [ ] **P1.5** 实现贪心细化（Greedy）
- [ ] **P1.6** 实现多层框架骨架（coarsen→partition→uncoarsen）
- [ ] **P1.7** 实现最简 2-way 划分通路（`PartGraphRecursive` k=2）
- [ ] **P1.8** 基础测试：正确性、平衡性断言

**完成标准**：`cargo test` 全通过，WASM 编译成功。

### Phase 2：算法完善

- [ ] **P2.1** 实现 SHEM 粗化
- [ ] **P2.2** 实现 FM 细化（2-way）
- [ ] **P2.3** 实现递归二分（任意 k）
- [ ] **P2.4** 实现 K-way 直接划分框架
- [ ] **P2.5** 实现 K-way FM 细化
- [ ] **P2.6** 实现多约束平衡（ncon > 1）
- [ ] **P2.7** 实现 `NodeND`（嵌套剖分 + 顶点分隔符）

**完成标准**：三个主函数完整可用，通过标准测试图的回归测试。

### Phase 3：API 层

- [ ] **P3.1** 实现 C ABI 兼容层（`c-api` feature）
- [ ] **P3.2** 实现 WASM/JS 绑定（`wasm` feature，wasm-bindgen）
- [ ] **P3.3** 编写 TypeScript 类型声明文件

### Phase 4：质量与性能

- [ ] **P4.1** 添加 fuzzing（cargo-fuzz）
- [ ] **P4.2** 添加 benchmark（criterion）
- [ ] **P4.3** 可选并行支持（rayon，`parallel` feature）
- [ ] **P4.4** WASM 尺寸优化（wasm-opt）

---

## 3. 模块实现指南

### 3.1 模块 `types.rs`

**职责**：定义所有共享类型。

```rust
// 实现要点：
// 1. Idx = i32，Real = f32（与 METIS 保持一致）
// 2. Options 实现 Default trait，对应 METIS_SetDefaultOptions 的默认值
// 3. MetisError 实现 thiserror::Error
// 4. Graph 实现 Display（调试用）

pub type Idx = i32;
pub type Real = f32;

// Options 默认值（来自 METIS 5.1.x 源码 options.c）：
// ctype   = SHEM
// iptype  = GROW（kway）/ METISRB（recursive）
// rtype   = GREEDY（kway）/ FM（recursive）
// ncuts   = 1
// niter   = 10
// seed    = -1（即随机）
// ufactor = 30（kway）/ 1（recursive）
// numbering = 0
```

### 3.2 模块 `graph/mod.rs`

**职责**：图的构造、验证、基本查询。

```rust
impl Graph {
    /// 构造函数：验证 CSR 格式合法性
    pub fn new(nvtxs: usize, xadj: Vec<Idx>, adjncy: Vec<Idx>, ...) -> Result<Self, MetisError>
    
    /// 顶点 v 的邻居切片
    pub fn neighbors(&self, v: usize) -> &[Idx]
    
    /// 顶点 v 的边权重切片（若存在）
    pub fn edge_weights(&self, v: usize) -> Option<&[Idx]>
    
    /// 顶点 v 的权重（多约束，返回 &[Idx] 长度 ncon）
    pub fn vertex_weights(&self, v: usize) -> &[Idx]
    
    /// 计算当前划分的边割数
    pub fn edge_cut(&self, part: &[Idx]) -> Idx
    
    /// 计算每个分区的权重向量
    pub fn partition_weights(&self, part: &[Idx], nparts: usize) -> Vec<Vec<Idx>>
    
    /// 提取边界顶点（邻居中有不同分区的顶点）
    pub fn boundary_vertices(&self, part: &[Idx]) -> Vec<usize>
}
```

**CSR 验证（validate.rs）必须检查：**
1. `xadj.len() == nvtxs + 1`
2. `xadj[0] == 0`，`xadj` 单调不减
3. `adjncy.len() == xadj[nvtxs]`
4. 所有 `adjncy[i]` 在 `[0, nvtxs)` 范围内
5. 无自环：`adjncy[j] != v` for all neighbors j of v
6. 对称性：`u` 是 `v` 的邻居则 `v` 也是 `u` 的邻居
7. 边权重数组长度与 adjncy 匹配

### 3.3 模块 `coarsen/`

**核心数据流：**

```
Graph  →  [Matcher]  →  matching: Vec<Idx>  →  [build_coarse_graph]  →  CoarseGraph
```

**关键实现细节：**

```rust
// coarsen/mod.rs
pub fn coarsen(graph: &Graph, options: &Options) -> Vec<CoarseGraph> {
    // 返回从细到粗的图层次（graph[0] 最细，graph.last() 最粗）
    let mut levels = vec![];
    let mut current = graph.clone();
    
    loop {
        let coarse = match options.ctype {
            CType::Rm   => rm::coarsen_step(&current),
            CType::Shem => shem::coarsen_step(&current),
        };
        
        // 停止条件
        let ratio = coarse.graph.nvtxs as f32 / current.nvtxs as f32;
        let threshold = (20 * nparts).max(options.coarsen_limit);
        
        levels.push(coarse.clone());
        
        if coarse.graph.nvtxs <= threshold || ratio > 0.95 {
            break;
        }
        current = coarse.graph;
    }
    levels
}

// build_coarse_graph 中的注意事项：
// 1. 超顶点权重 = 两个匹配顶点权重之和（多约束时每个维度独立相加）
// 2. 边去重时：同一对超顶点间的多条边合并为一条，权重相加
// 3. 保存 cmap：原始顶点 v → 粗化超顶点 ID
// 4. 不匹配的顶点（孤立顶点）自映射到单独超顶点
```

### 3.4 模块 `initial/`

**grow.rs 实现要点：**

```rust
// 图增长法（BFS-based）
pub fn grow_partition(graph: &Graph, nparts: usize, tpwgts: &[Real], rng: &mut impl Rng) -> Vec<Idx> {
    // 1. 随机选择 nparts 个种子顶点
    // 2. 对每个种子做 BFS，直到达到目标权重
    // 3. 关键：随机化 BFS 队列处理顺序，避免形状偏差
    // 4. 剩余未分配顶点贪心分配给权重最小的分区
    
    // 目标权重：tpwgts[p] × total_weight（对每个约束独立计算）
    // 分配优先级：边权重加权的邻居优先（SHEM 精神）
}
```

### 3.5 模块 `refine/`

**fm.rs（2-way FM）实现要点：**

```rust
pub struct FmRefiner<'a> {
    graph: &'a Graph,
    part: Vec<Idx>,
    // 每个顶点的增益：移动到另一侧可减少的边割数
    gain: Vec<Idx>,
    // 边界顶点集合
    boundary: HashSet<usize>,
    // 每个分区的当前权重
    part_weight: [Idx; 2],
}

impl<'a> FmRefiner<'a> {
    // 计算顶点 v 的初始增益
    fn compute_gain(&self, v: usize) -> Idx {
        // gain = (edges to other part) - (edges to same part)
        let mut gain = 0;
        for (u, w) in self.graph.weighted_neighbors(v) {
            if self.part[u] == self.part[v] {
                gain -= w;  // 内部边，移动后变为切割边
            } else {
                gain += w;  // 切割边，移动后变为内部边
            }
        }
        gain
    }
    
    // 执行一次 FM pass
    // 返回：是否有改进
    fn fm_pass(&mut self, niter: usize, ubvec: f32) -> bool {
        // 使用 priority-queue crate 的 PriorityQueue<usize, Idx>
        // 关键：移动后更新所有邻居的增益
        // 回滚机制：记录所有移动，回滚到最佳切割点
    }
}
```

**fm_kway.rs（K-way FM）扩展：**

```rust
// 每个顶点对每个分区维护增益
// gain[v][p] = 移动 v 到分区 p 的增益（相对于当前分区）
// 内存优化：只对边界顶点维护增益，非边界顶点增益为 -∞
// 数据结构：Vec<HashMap<PartId, Idx>> 或 Vec<BTreeMap<PartId, Idx>>
```

### 3.6 模块 `partition/`

**recursive.rs：**

```rust
pub fn partition_recursive(
    graph: &Graph,
    nparts: usize,
    tpwgts: &[Real],
    ubvec: &[Real],
    options: &Options,
    rng: &mut impl Rng,
) -> Result<PartitionResult, MetisError> {
    // 递归二分实现：
    // 1. 基础情况：nparts == 1，全部分到 0 号分区
    // 2. 分割：k1 = nparts/2，k2 = nparts - k1
    // 3. 对完整多层框架执行二分
    // 4. 递归处理两个子图（注意：子图的顶点 ID 需要重映射）
    // 5. 合并结果时偏移第二半的分区 ID（+= k1）
    
    // 重要：子图构造时必须重建 CSR，保持紧凑顶点编号
}
```

**kway.rs：**

```rust
pub fn partition_kway(
    graph: &Graph,
    nparts: usize,
    ...
) -> Result<PartitionResult, MetisError> {
    // 1. 多层粗化（同 recursive）
    // 2. 初始划分：直接在粗化图上做 k-way 初始划分
    //    - 使用 grow_partition 或 random_partition
    // 3. 细化：使用 K-way FM 或 Greedy
    // 4. 投影：将粗化图划分逐层投影回原始图
    
    // 投影算法：
    // for each level from coarsest to finest:
    //   for each vertex v in fine graph:
    //     part[v] = coarse_part[cmap[v]]
    //   refine on fine graph
}
```

### 3.7 模块 `separator/mod.rs`

**用于 NodeND 的顶点分隔符：**

```rust
// 算法：基于边分隔符转换为顶点分隔符
// 步骤：
// 1. 对图做二分（使用 multilevel bisection）
// 2. 找到跨越两个分区的边
// 3. 从较小分区选择顶点组成分隔符
// 4. 确保移除分隔符后图不再连通

pub fn find_node_separator(
    graph: &Graph,
    options: &Options,
    rng: &mut impl Rng,
) -> (Vec<usize>, Vec<usize>, Vec<usize>)
// 返回：(left_part, right_part, separator)
```

### 3.8 模块 `ffi/c_api.rs`

```rust
// 注意：这是唯一允许使用 unsafe 的模块
// 所有 C ABI 函数必须：
// 1. 在函数入口验证所有指针非空
// 2. 将 Rust 错误转换为 METIS 错误码
// 3. 不 panic（使用 catch_unwind 包装或避免 panic 路径）

#[no_mangle]
pub unsafe extern "C" fn RMETIS_PartGraphKway(
    nvtxs: *const Idx, ncon: *const Idx,
    xadj: *const Idx, adjncy: *const Idx,
    vwgt: *const Idx, vsize: *const Idx, adjwgt: *const Idx,
    nparts: *const Idx, tpwgts: *const Real, ubvec: *const Real,
    options: *const Idx,
    objval: *mut Idx, part: *mut Idx,
) -> i32 {
    // 1. 空指针检查（必须的 nvtxs, ncon, xadj, adjncy, nparts, part）
    // 2. 从裸指针构建 Rust 切片（unsafe）
    // 3. 调用安全的 Rust API
    // 4. 写入结果到 objval 和 part
    // 5. 错误映射到 METIS_ERROR_* 常量
}
```

### 3.9 模块 `ffi/wasm_api.rs`

```rust
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct WasmGraph { /* ... */ }

#[wasm_bindgen]
impl WasmGraph {
    #[wasm_bindgen(constructor)]
    pub fn new(nvtxs: usize, xadj: &[i32], adjncy: &[i32]) -> Result<WasmGraph, JsValue>
    
    pub fn set_vertex_weights(&mut self, vwgt: &[i32], ncon: usize)
    pub fn set_edge_weights(&mut self, adjwgt: &[i32])
}

#[wasm_bindgen]
pub fn part_graph_kway(
    graph: &WasmGraph,
    nparts: usize,
    options_json: Option<String>,  // JSON 字符串传递 Options
) -> Result<JsValue, JsValue>  // 返回 {part: Int32Array, objval: number}
```

---

## 4. 关键算法实现注意事项

### 4.1 优先队列选择

FM 算法依赖高效的优先队列（需要支持按 key 查找 + 修改优先级）：

```toml
# 推荐使用 priority-queue crate
priority-queue = "1.3"
```

```rust
use priority_queue::PriorityQueue;
// PriorityQueue<VertexId, Gain>
// 支持 O(log n) 的 push, pop_max, change_priority
let mut pq: PriorityQueue<usize, i32> = PriorityQueue::new();
pq.push(vertex_id, gain);
pq.change_priority(&vertex_id, new_gain);
let (best_vertex, best_gain) = pq.pop().unwrap();
```

### 4.2 随机数生成

```rust
use rand::{SeedableRng, Rng};
use rand::rngs::SmallRng;

// 从 Options.seed 创建 RNG
let rng = if options.seed == 0 {
    // WASM 下用 getrandom，native 下用 thread_rng seed
    SmallRng::from_entropy()
} else {
    SmallRng::seed_from_u64(options.seed as u64)
};
```

### 4.3 整数溢出防范

METIS 使用 `i32`，大图的权重累加可能溢出：

```rust
// 使用饱和加法或检查溢出
let total_weight: i64 = vwgt.iter().map(|&w| w as i64).sum();
assert!(total_weight <= i32::MAX as i64, "vertex weight sum overflow");
```

### 4.4 多约束权重布局

```
// 内存布局：行优先，顶点为行，约束为列
// vwgt[v * ncon + c] = 顶点 v 的第 c 个约束权重
// tpwgts[p * ncon + c] = 分区 p 的第 c 个约束目标权重

// 获取顶点 v 的权重向量：
fn vertex_weight_slice(vwgt: &[Idx], v: usize, ncon: usize) -> &[Idx] {
    &vwgt[v * ncon .. (v + 1) * ncon]
}
```

### 4.5 子图提取（递归二分用）

```rust
// 从划分结果提取子图，重建紧凑 CSR
fn extract_subgraph(graph: &Graph, part: &[Idx], target_part: Idx) -> (Graph, Vec<usize>) {
    // 返回：(子图, 子图顶点到原始图顶点的映射)
    let vertices: Vec<usize> = (0..graph.nvtxs)
        .filter(|&v| part[v] == target_part)
        .collect();
    
    let old_to_new: HashMap<usize, usize> = vertices.iter()
        .enumerate()
        .map(|(new, &old)| (old, new))
        .collect();
    
    // 重建 CSR（只保留子图内部的边）
    // ...
}
```

### 4.6 投影（Projection）

```rust
// 从粗化图划分投影到细化图
fn project_partition(
    fine_graph: &Graph,
    coarse_part: &[Idx],
    cmap: &[Idx],  // fine顶点 → coarse超顶点
) -> Vec<Idx> {
    (0..fine_graph.nvtxs)
        .map(|v| coarse_part[cmap[v] as usize])
        .collect()
}
```

---

## 5. 测试策略

### 5.1 单元测试结构

每个算法模块在 `#[cfg(test)]` 块中包含：

```rust
#[cfg(test)]
mod tests {
    use super::*;
    
    fn small_graph() -> Graph {
        // 构造标准测试图（见附录 B）
    }
    
    #[test]
    fn test_csr_validity() { ... }
    
    #[test]
    fn test_partition_balance() {
        // 验证每个分区权重 ≤ (1 + ufactor/1000) × total/nparts
    }
    
    #[test]
    fn test_partition_coverage() {
        // 验证所有顶点被分配（part[v] ∈ [0, nparts)）
    }
}
```

### 5.2 集成测试（tests/correctness.rs）

```rust
// 使用标准 METIS 测试图
// 测试目录：tests/fixtures/
// 文件格式：METIS graph format（.graph 文件）

#[test]
fn test_4elt_kway_k8() {
    let graph = load_metis_graph("tests/fixtures/4elt.graph");
    let result = part_graph_kway(&graph, 8, None, None, &Options::default()).unwrap();
    
    // 验证正确性
    assert_partition_valid(&graph, &result, 8);
    
    // 验证平衡性（容差 30/1000 = 3%）
    assert_partition_balanced(&graph, &result, 8, 30);
}
```

### 5.3 回归测试（tests/regression.rs）

```rust
// 固定 seed，验证边割数的确定性
#[test]
fn test_deterministic_output() {
    let graph = load_metis_graph("tests/fixtures/4elt.graph");
    let opts = Options { seed: 42, ..Default::default() };
    
    let r1 = part_graph_kway(&graph, 4, None, None, &opts).unwrap();
    let r2 = part_graph_kway(&graph, 4, None, None, &opts).unwrap();
    
    assert_eq!(r1.objval, r2.objval);
    assert_eq!(r1.part, r2.part);
}
```

### 5.4 辅助断言函数

```rust
// tests/common.rs

pub fn assert_partition_valid(graph: &Graph, result: &PartitionResult, nparts: usize) {
    assert_eq!(result.part.len(), graph.nvtxs);
    for &p in &result.part {
        assert!(p >= 0 && (p as usize) < nparts, "invalid partition id: {}", p);
    }
}

pub fn assert_partition_balanced(
    graph: &Graph, result: &PartitionResult, nparts: usize, ufactor: u32
) {
    let weights = graph.partition_weights(&result.part, nparts);
    let total: Idx = graph.total_weight();
    let target = total as f32 / nparts as f32;
    let max_allowed = target * (1.0 + ufactor as f32 / 1000.0);
    
    for (p, pw) in weights.iter().enumerate() {
        let w = pw[0] as f32;
        assert!(w <= max_allowed,
            "partition {} weight {} exceeds max allowed {}", p, w, max_allowed);
    }
}

pub fn assert_nd_valid(nvtxs: usize, result: &NDResult) {
    // 验证 perm 和 iperm 互为逆置换
    let mut check = vec![false; nvtxs];
    for &p in &result.perm {
        assert!(!check[p as usize], "duplicate in perm");
        check[p as usize] = true;
    }
    for v in 0..nvtxs {
        assert_eq!(result.perm[result.iperm[v] as usize] as usize, v);
    }
}
```

---

## 6. WASM 构建流程

### 6.1 依赖工具

```bash
# 安装 wasm-pack
cargo install wasm-pack

# 安装 wasm-opt（Binaryen）
# macOS: brew install binaryen
# Linux: apt install binaryen
# Windows: 从 https://github.com/WebAssembly/binaryen/releases 下载
```

### 6.2 构建命令

```bash
# 开发构建（调试）
wasm-pack build --target web --dev -- --features wasm

# 生产构建
wasm-pack build --target web --release -- --features wasm

# Node.js 目标
wasm-pack build --target nodejs --release -- --features wasm

# 手动 wasm-opt 优化（wasm-pack 会自动调用，但可以手动指定级别）
wasm-opt -O3 -o pkg/rmetis_bg_opt.wasm pkg/rmetis_bg.wasm
```

### 6.3 WASM 测试

```bash
# 在无头浏览器中运行测试
wasm-pack test --headless --firefox -- --features wasm

# Node.js 环境测试
wasm-pack test --node -- --features wasm
```

### 6.4 WASM 使用示例

```html
<!-- index.html -->
<script type="module">
import init, { WasmGraph, part_graph_kway } from './pkg/rmetis.js';

await init();

// 构建图：0-1, 1-2, 2-3, 3-0 (正方形)
const graph = new WasmGraph(4,
  new Int32Array([0, 2, 4, 6, 8]),   // xadj
  new Int32Array([1, 3, 0, 2, 1, 3, 0, 2])  // adjncy
);

const result = part_graph_kway(graph, 2);
console.log('partition:', result.part);      // Int32Array
console.log('edge cut:', result.objval);     // number
</script>
```

---

## 7. C ABI 构建与使用

### 7.1 构建共享库

```toml
# Cargo.toml
[lib]
crate-type = ["cdylib", "rlib"]
```

```bash
# 构建动态库
cargo build --release --features c-api

# Linux: target/release/librmetis.so
# macOS: target/release/librmetis.dylib
# Windows: target/release/rmetis.dll
```

### 7.2 头文件生成

提供 `include/rmetis.h`：

```c
#ifndef RMETIS_H
#define RMETIS_H

#include <stdint.h>

typedef int32_t idx_t;
typedef float   real_t;

#define RMETIS_OK           1
#define RMETIS_ERROR_INPUT -2
#define RMETIS_ERROR_MEMORY -3
#define RMETIS_ERROR       -4

#define RMETIS_NOPTIONS    40

int RMETIS_SetDefaultOptions(idx_t *options);
int RMETIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy,
    idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts,
    real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, idx_t *part);
int RMETIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy,
    idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts,
    real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, idx_t *part);
int RMETIS_NodeND(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
    idx_t *options, idx_t *perm, idx_t *iperm);

#endif
```

---

## 8. CI/CD 配置

### 8.1 GitHub Actions 工作流

```yaml
# .github/workflows/ci.yml
name: CI

on: [push, pull_request]

jobs:
  test:
    strategy:
      matrix:
        os: [ubuntu-latest, windows-latest, macos-latest]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo test --all-features
      - run: cargo clippy -- -D warnings
      - run: cargo fmt --check

  wasm:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          targets: wasm32-unknown-unknown
      - run: cargo install wasm-pack
      - run: wasm-pack build --target web --release -- --features wasm
      - run: wasm-pack test --node -- --features wasm

  bench:
    runs-on: ubuntu-latest
    if: github.event_name == 'push' && github.ref == 'refs/heads/main'
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo bench --no-run  # 仅编译不运行（PR 阶段）
```

---

## 9. 文件格式规范

### 9.1 METIS Graph 格式（.graph 文件）

用于测试数据输入：

```
% 注释行以 % 开头
nvtxs nedges [fmt] [ncon]
# fmt: 0=无权重，1=有边权重，10=有顶点权重，11=均有
# 以下每行为一个顶点的邻居列表（空格分隔）
# fmt=1: 邻居 权重 邻居 权重 ...
# fmt=10: 顶点权重 邻居 邻居 ...
# fmt=11: 顶点权重 邻居 权重 邻居 权重 ...

# 示例：4顶点，4边，无权重
4 4
2 4
1 3
2 4
1 3
```

### 9.2 图格式读取器（tests/common.rs）

```rust
pub fn load_metis_graph(path: &str) -> Graph {
    // 解析 METIS .graph 格式
    // 注意：METIS 文件使用 1-based 编号，转换为 0-based
}
```

---

## 10. 性能分析指导

### 10.1 热点函数

预期性能瓶颈（按顺序）：

1. `fm_pass`：FM 细化的优先队列操作
2. `shem_coarsen`：邻居边权重排序
3. `build_coarse_graph`：边去重（HashMap 操作）
4. `compute_initial_gains`：FM 初始化

### 10.2 性能优化方向（v0.2）

- 用 `IndexSet` 替代 `HashSet` 提升边界顶点遍历缓存局部性
- 粗化图构造用 `Vec<(Idx, Idx)>` + sort_unstable 代替 HashMap 去重
- FM 增益数组预分配，避免频繁 realloc
- K-way FM 使用 `SmallVec` 存储每个顶点的邻接分区增益

### 10.3 Benchmark 结构

```rust
// benches/partition.rs
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

fn bench_kway(c: &mut Criterion) {
    let graphs = [
        ("4elt", load_metis_graph("tests/fixtures/4elt.graph")),
        ("mdual", load_metis_graph("tests/fixtures/mdual.graph")),
    ];
    
    let mut group = c.benchmark_group("PartGraphKway");
    for (name, graph) in &graphs {
        for nparts in [4, 8, 16, 32] {
            group.bench_with_input(
                BenchmarkId::new(*name, nparts),
                &(graph, nparts),
                |b, (g, k)| b.iter(|| {
                    part_graph_kway(g, *k, None, None, &Options::default())
                }),
            );
        }
    }
}

criterion_group!(benches, bench_kway);
criterion_main!(benches);
```

---

## 附录 A：开发顺序建议（AI Agent 执行顺序）

AI Agent 按以下顺序依次实现，每个步骤完成后运行 `cargo test` 确认无回归：

```
1. types.rs          → 类型和错误定义
2. graph/validate.rs → CSR 验证逻辑
3. graph/mod.rs      → Graph 结构体和基本方法
4. coarsen/rm.rs     → 随机匹配粗化
5. initial/grow.rs   → 图增长初始划分
6. refine/greedy.rs  → 贪心细化
7. partition/kway.rs → K-way 多层框架（先用贪心细化）
8. tests/correctness → 第一批集成测试
9. coarsen/shem.rs   → SHEM 粗化（替换 RM）
10. refine/fm.rs     → 2-way FM 细化
11. partition/recursive.rs → 递归二分
12. refine/fm_kway.rs → K-way FM
13. separator/mod.rs  → 顶点分隔符
14. partition/nd.rs   → NodeND
15. ffi/c_api.rs      → C ABI 层
16. ffi/wasm_api.rs   → WASM 绑定
```

## 附录 B：标准测试图

**小型验证图（手工构造，写入代码）：**

```rust
// 4x4 网格图（16顶点，24边）
fn grid4x4() -> Graph {
    // 顶点编号：
    // 0  1  2  3
    // 4  5  6  7
    // 8  9  10 11
    // 12 13 14 15
    // 边：水平 + 垂直相邻
}

// 路径图（10顶点）
fn path10() -> Graph { ... }

// 完全二叉树（15顶点，深度4）
fn binary_tree() -> Graph { ... }
```

**标准数据集（从 METIS 官方获取）：**

- `tests/fixtures/4elt.graph`：15,606 顶点，45,878 边（有限元网格）
- `tests/fixtures/copter2.graph`：9,167 顶点，49,152 边（直升机结构）

## 附录 C：常见陷阱

| 陷阱 | 描述 | 解决方案 |
|------|------|---------|
| METIS 1-based | .graph 文件用 1-based 编号 | 读取时统一减 1 转 0-based |
| 无向图双向存储 | CSR 中每条无向边存两次 | nedges = adjncy.len() / 2 |
| tpwgts 归一化 | 必须对每个 ncon 独立归一化（和为1） | 加载时验证或归一化 |
| FM 增益更新 | 移动顶点后必须更新所有邻居增益 | 不要忘记更新锁定顶点的邻居 |
| 粗化图节点数 | 孤立匹配顶点也必须成为超顶点 | matching[-1] 的顶点自映射 |
| WASM 无 panic | wasm32 的 panic 会终止整个 wasm 实例 | 所有公开函数用 Result 返回，不 unwrap |
