import React, { useState } from 'react';
import { analyzeFile } from '../api';

function UploadForm({ setResult }) {
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!file) return;
    setLoading(true);
    const res = await analyzeFile(file);
    setResult(res);
    setLoading(false);
  };

  return (
    <form onSubmit={handleSubmit}>
      <input type="file" accept=".h5ad,.csv" onChange={e => setFile(e.target.files[0])} />
      <button type="submit" disabled={loading}>{loading ? '分析中...' : '上传并分析'}</button>
    </form>
  );
}

export default UploadForm;