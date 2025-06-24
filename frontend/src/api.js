export async function analyzeFile(file) {
  const formData = new FormData();
  formData.append('file', file);
  const res = await fetch('http://localhost:8000/analyze/', {
    method: 'POST',
    body: formData,
  });
  return await res.json();
}