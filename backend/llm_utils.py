import json
import logging
import subprocess
import uuid
import os
from scanpy_utils import plot_violin, plot_qc, plot_umap, check_adata_structure, plot_dotplot, plot_pseudotime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ANALYSIS_FUNCTIONS = {
    "plot_violin": plot_violin,
    "plot_qc": plot_qc,
    "plot_umap": plot_umap,
    "check_adata_structure": check_adata_structure,
    "plot_dotplot": plot_dotplot,
    "plot_pseudotime": plot_pseudotime
}

def get_prompt(user_request, adata=None):
    genes = list(adata.var_names[:10]) if adata is not None else []
    
    # 直接使用固定的分组变量列表
    allowed_groups = ["group", "celltype", "orig.ident"]
    
    return (
    f"你是单细胞分析助手。解析用户请求并输出JSON格式的命令。\n"
    f"支持的功能有：\n"
    f"1. plot_violin：绘制基因表达的小提琴图\n"
    f"2. plot_umap：绘制UMAP降维可视化\n" 
    f"3. plot_qc：绘制QC质量控制图\n"
    f"4. check_adata_structure：检查数据结构\n"
    f"5. plot_dotplot：绘制点图\n"
    f"6. plot_pseudotime：绘制拟时序分析\n\n"
    
    f"重要！必须严格遵守以下分组规则，不允许任何例外：\n"
    f"- 无论用户如何描述，只要是关于比较CP和Ctrl或任何分组比较的请求，都必须使用groupby=\"group\"\n"
    f"- 当用户提到按照细胞类型绘图时，必须使用groupby=\"celltype\"\n"
    f"- 当用户提到按照细胞簇绘图时，必须使用groupby=\"orig.ident\"\n"
    f"- 永远不要使用nCount_RNA, nFeature_RNA等其他列名\n"
    f"- 不要创建新的分组变量\n"
    f"可用基因示例: {genes},不区分大小写，收到基因名后统一改成第一个字母大写\n"
    f"唯一可用的分组变量: group, celltype, orig.ident\n"
    f"用户输入：{user_request}\n\n"
    
    f"返回格式：\n"
    f"{{\"function\": \"plot_violin\", \"params\": {{\"genes\": [\"基因1\", \"基因2\"], \"groupby\": \"分组变量\"}}}}\n"
    f"或：{{\"function\": \"check_adata_structure\", \"params\": {{}}}}\n"
)

def extract_first_json_block(text):
    start = text.find('{')
    if start == -1:
        return None
    count = 0
    for i in range(start, len(text)):
        if text[i] == '{':
            count += 1
        elif text[i] == '}':
            count -= 1
        if count == 0:
            return text[start:i+1]
    return None

class RemoteOllamaLLM:
    def __init__(self, model="qwen3:4b"):
        self.model = model

    def __call__(self, user_request, adata=None):
        prompt = get_prompt(user_request, adata)
        try:
            cmd = ["ollama", "run", self.model, prompt]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            content = result.stdout
            logger.info(f"LLM raw response: {content}")

            snippet = extract_first_json_block(content)
            if not snippet:
                logger.error("No JSON block found in LLM response")
                return {"error": "no json found in LLM response"}
            logger.info(f"Extracted JSON snippet: {snippet!r}")

            try:
                parsed = json.loads(snippet)
                if isinstance(parsed, dict) and "function" in parsed:
                    logger.info(f"Extracted valid JSON: {snippet!r}")
                    return {"function": parsed.get("function"), "params": parsed.get("params", {})}
                else:
                    logger.error("No function field in parsed JSON")
                    return {"error": "no function field in LLM response"}
            except Exception as e:
                logger.error(f"JSON decode error: {e}, snippet: {snippet!r}")
                return {"error": "invalid json in LLM response,shit man"}
        except Exception as e:
            logger.error(f"LLM or JSON error: {e}")
            return {"error": f"LLM/JSON error: {str(e)}"}

import uuid
def parse_and_execute(user_request, llm_model, adata):
    parsed_output = llm_model(user_request, adata)
    logger.info(f"Parsed LLM output: {parsed_output}")
    if "error" in parsed_output:
        return parsed_output

    function_name = parsed_output.get("function")
    params = parsed_output.get("params", {})

    # 参数名标准化映射
    param_aliases = {
        "genes": ["genes", "gene", "gene_list", "features"],
        "groupby": ["groupby", "group", "group_by"],
    }
    for std_key, aliases in param_aliases.items():
        for alias in aliases:
            if alias in params:
                if std_key == "genes" and isinstance(params[alias], str):
                    params[std_key] = [params[alias]]
                else:
                    params[std_key] = params[alias]
    
        # 兜底：如果groupby/group/group_by为seurat相关，直接移除或替换为默认分组
    seurat_exclude = {"seurat_clusters", "orig.ident", "group"}
    # 参数名标准化映射
    param_aliases = {
        "genes": ["genes", "gene", "gene_list", "features"],
        "groupby": ["groupby", "group", "group_by"],
    }
    for std_key, aliases in param_aliases.items():
        for alias in aliases:
            if alias in params:
                # genes 相关参数转为 list
                if std_key == "genes" and isinstance(params[alias], str):
                    params[std_key] = [params[alias]]
                else:
                    params[std_key] = params[alias]

    # groupby 只允许单个分组且不能是 seurat 相关
    seurat_exclude = {"seurat_clusters", "orig.ident", "group"}
    if "groupby" in params:
        # 如果是列表，取第一个合法分组
        if isinstance(params["groupby"], list):
            candidates = [g for g in params["groupby"] if g not in seurat_exclude and g in adata.obs.columns]
            params["groupby"] = candidates[0] if candidates else None
        # 如果是非法分组，自动换成第一个合法分组
        elif params["groupby"] in seurat_exclude or params["groupby"] not in adata.obs.columns:
            candidates = [g for g in adata.obs.columns if g not in seurat_exclude]
            params["groupby"] = candidates[0] if candidates else adata.obs.columns[0]
    # 如果没有groupby，自动补一个合法分组
    if "groupby" not in params or not params["groupby"]:
        candidates = [g for g in adata.obs.columns if g not in seurat_exclude]
        params["groupby"] = candidates[0] if candidates else adata.obs.columns[0]
    func = ANALYSIS_FUNCTIONS[function_name]
    import inspect
    sig = inspect.signature(func)
    valid_keys = set(sig.parameters.keys())
    filtered_params = {k: v for k, v in params.items() if k in valid_keys}

    base_path = f"static/{uuid.uuid4().hex}"
    pdf_path = base_path + ".pdf"
    png_path = base_path + ".png"
    tiff_path = base_path + ".tiff"

    # 生成pdf
    filtered_params["output_path"] = pdf_path
    func(adata, **filtered_params)

    # 生成png
    from scanpy_utils import plot_violin, plot_qc, plot_umap, check_adata_structure, plot_dotplot, plot_pseudotime
    filtered_params["output_path"] = png_path
    func(adata, **filtered_params)

    # 生成tiff
    filtered_params["output_path"] = tiff_path
    func(adata, **filtered_params)

    return {
        "pdf_url": "/" + pdf_path,
        "png_url": "/" + png_path,
        "tiff_url": "/" + tiff_path
    }