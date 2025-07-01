from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
import scanpy as sc
import tempfile
import os
from llm_utils import parse_and_execute, RemoteOllamaLLM
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")
app.mount("/frontend", StaticFiles(directory="../frontend/src"), name="frontend")

# 静态文件目录
STATIC_DIR = "static"
os.makedirs(STATIC_DIR, exist_ok=True)


# 跨域
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 初始化 LLM
llm_model = RemoteOllamaLLM(model="qwen3:4b")

@app.get("/")
def read_index():
    return FileResponse("../frontend/src/index.html")

@app.get("/hello_qwen/")
async def hello_qwen():
    try:
        return llm_model("hello")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/analyze/")
async def analyze(file: UploadFile = File(...), user_request: str = Form(...)):
    tmp_path = None
    try:
        content = await file.read()
        print("Received file size:", len(content), "bytes")
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(content)
            tmp_path = tmp.name

        # 读取 AnnData
        adata = sc.read(tmp_path)

        # 处理请求
        result = parse_and_execute(user_request, llm_model, adata)

        # 返回图片链接或结果
        if isinstance(result, dict) and "image_url" in result:
            return result
        return result

    except Exception as e:
        print("Error:", str(e))
        return {"error": str(e)}
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)