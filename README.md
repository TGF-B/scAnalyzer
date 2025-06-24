# scAnalyzer

A modern, SpaceX-style web application for single-cell data analysis powered by Qwen-7B.

## Features

- **Upload AnnData (.h5ad) files** for analysis
- **Flexible user requests**: supports violin plots, UMAP, dotplot, QC, pseudotime, and more
- **Interactive web frontend** with stylish UI
- **Download results** in PDF, PNG, and TIFF formats
- **Backend powered by FastAPI and Scanpy**
- **LLM-driven analysis**: leverages Qwen-7B for intelligent request parsing

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/TGF-B/scAnalyzer.git
cd scAnalyzer
```

### 2. Install backend dependencies

```bash
cd backend
pip install -r requirements.txt
```

### 3. Start the backend server

```bash
uvicorn app:app --host 0.0.0.0 --port 8005
```

### 4. Serve the frontend

You can use any static file server, or integrate with your backend as needed.

### 5. Open the app

Visit [http://localhost:8005](http://localhost:8005) in your browser.

## Directory Structure

```
scAnalyzer/
├── backend/         # FastAPI backend and analysis logic
├── frontend/        # Frontend (HTML/CSS/JS, React or static)
├── .gitignore
└── README.md
```

## .gitignore Highlights

- Ignores Python cache, large data files, images, and node_modules
- Keeps your repository clean and fast

## Notes

- **Do not commit large data files (.h5ad, images, etc.)** — keep them out of git for performance.
- For deployment, make sure to configure your static file serving and backend endpoints as needed.

## License

MIT

---

**Thanks for using scAnalyzer! If you have questions or suggestions, feel free to open an issue or pull request.**
