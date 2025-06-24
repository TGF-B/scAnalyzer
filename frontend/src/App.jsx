import React, { useState } from 'react';
import UploadForm from './components/UploadForm';
import ResultDisplay from './components/ResultDisplay';

function App() {
  const [result, setResult] = useState(null);

  return (
    <div style={{ padding: 32 }}>
      <h1>单细胞分析大模型平台</h1>
      <UploadForm setResult={setResult} />
      <ResultDisplay result={result} />
    </div>
  );
}

export default App;