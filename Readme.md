# GenoCOPoint

## 1. 简介
GenoCOPoint 是一个用于检测基因型数据中变点（CO点）的R包。它支持读取五列格式的基因型文件，结合真实CO区间文件，自动检测变点并生成可视化PDF和结果。

## 2. 依赖
- R (>= 3.5.0)
- [changepoint](https://github.com/rkillick/changepoint/) R包

## 3. 安装
- 详情请见[Installation.md](https://github.com/StellarMech/GenoCOPoint/blob/main/Installation.md)

## 4. 使用方法

### 4.1 函数说明

主函数为 `GenoCOPoint()`，参数如下：

- `geno_file`：基因型文件路径（五列格式）
- `truth_file`：真实CO区间文件路径
- `output_pdf`：输出PDF文件名
- `output_txt`：输出结果表格文件名
- `window_size`：平滑窗口大小
- `smooth_threshold`：平滑阈值
- `min_support`：最小连续SNP数
- `consistency`：一致性比例

### 4.2 示例

```r
library(GenoCOPoint)
geno_file <- "GenoCOPoint/examples/genofiles/AAACGCTCATATAGCC_smt_genotypes_withname_reads.txt"
truth_file <- "GenoCOPoint/examples/truthfiles/AAACGCTCATATAGCC_allele_cnts_at_markers_sorted_co_pred.txt"
result <- GenoCOPoint(
  geno_file,
  truth_file,
  output_pdf = "example_output.pdf",
  output_txt = "example_output.txt",
  window_size = 5,
  smooth_threshold = 0.6,
  min_support = 7,
  consistency = 0.85
)
print(result)
```

## 5. 结果说明

- 输出PDF文件：每条染色体若有检测点两图展示（Reads柱状图+基因型散点图，叠加真实CO区间），无检测点一图展示（Reads柱状图）
- 输出TXT文件：每条染色体检测到的CO点位置信息

---

如需更多帮助，请参考包内文档或联系作者。