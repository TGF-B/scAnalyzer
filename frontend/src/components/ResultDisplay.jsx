import React from 'react';

function ResultDisplay({ result }) {
  if (!result) return null;
  return (
    <div style={{ marginTop: 24 }}>
      <h2>分析结果</h2>
      <pre>{JSON.stringify(result, null, 2)}</pre>
    </div>
  );
}

export default ResultDisplay;