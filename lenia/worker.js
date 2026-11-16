// 3D Lenia core: FFT convolution, single channel, multi-ring kernel
function makeFFT(n){
  const levels = Math.log2(n)|0; if ((1<<levels)!==n) throw 'pow2';
  const rev = new Uint32Array(n);
  for (let i=0;i<n;i++){ let r=0; for(let b=0;b<levels;b++) r|=((i>>b)&1)<<(levels-1-b); rev[i]=r; }
  const cos=new Float64Array(n/2), sin=new Float64Array(n/2);
  for(let i=0;i<n/2;i++){cos[i]=Math.cos(2*Math.PI*i/n); sin[i]=Math.sin(2*Math.PI*i/n);}
  function fft(re, im, inv){
    for(let i=0;i<n;i++){const j=rev[i]; if(j>i){let t=re[i];re[i]=re[j];re[j]=t;t=im[i];im[i]=im[j];im[j]=t;}}
    for(let size=2;size<=n;size<<=1){const half=size>>1, step=n/size;
      for(let i=0;i<n;i+=size){ for(let j=0,k=0;j<half;j++,k+=step){
        const c=cos[k], s=inv?sin[k]:-sin[k];
        const a=i+j, b=a+half;
        const tre=re[b]*c-im[b]*s, tim=re[b]*s+im[b]*c;
        re[b]=re[a]-tre; im[b]=im[a]-tim; re[a]+=tre; im[a]+=tim; }}}
  }
  return fft;
}
function makeLenia(N){
  const NC=N*N*N, N2=N*N;
  const fft=makeFFT(N);
  const re=new Float32Array(NC), im=new Float32Array(NC);
  const lr=new Float64Array(N), li=new Float64Array(N);
  let Kf=new Float32Array(NC);
  function fft3(inv){
    // x lines
    for(let z=0;z<N;z++)for(let y=0;y<N;y++){const b=z*N2+y*N;
      for(let x=0;x<N;x++){lr[x]=re[b+x];li[x]=im[b+x];} fft(lr,li,inv);
      for(let x=0;x<N;x++){re[b+x]=lr[x];im[b+x]=li[x];}}
    for(let z=0;z<N;z++)for(let x=0;x<N;x++){const b=z*N2+x;
      for(let y=0;y<N;y++){lr[y]=re[b+y*N];li[y]=im[b+y*N];} fft(lr,li,inv);
      for(let y=0;y<N;y++){re[b+y*N]=lr[y];im[b+y*N]=li[y];}}
    for(let y=0;y<N;y++)for(let x=0;x<N;x++){const b=y*N+x;
      for(let z=0;z<N;z++){lr[z]=re[b+z*N2];li[z]=im[b+z*N2];} fft(lr,li,inv);
      for(let z=0;z<N;z++){re[b+z*N2]=lr[z];im[b+z*N2]=li[z];}}
  }
  const core = r => (r<=0||r>=1)?0:Math.exp(4-1/(r*(1-r)));
  function setKernel(R, beta){
    re.fill(0); im.fill(0); let sum=0; const B=beta.length;
    for(let z=0;z<N;z++)for(let y=0;y<N;y++)for(let x=0;x<N;x++){
      const dx=x<N/2?x:x-N, dy=y<N/2?y:y-N, dz=z<N/2?z:z-N;
      const r=Math.sqrt(dx*dx+dy*dy+dz*dz)/R; if(r>=1) continue;
      const br=B*r, k=Math.min(B-1,Math.floor(br)); const v=beta[k]*core(br-k);
      re[z*N2+y*N+x]=v; sum+=v; }
    for(let i=0;i<NC;i++) re[i]/=sum;
    fft3(false); Kf=re.slice(); // kernel symmetric -> real spectrum
  }
  // A: Float32Array state; U written into U
  function step(A, U, mu, sigma, dt){
    re.set(A); im.fill(0); fft3(false);
    for(let i=0;i<NC;i++){re[i]*=Kf[i]; im[i]*=Kf[i];}
    fft3(true);
    const inv=1/NC, s2=2*sigma*sigma; let mass=0;
    for(let i=0;i<NC;i++){ const u=re[i]*inv; U[i]=u; const g=2*Math.exp(-(u-mu)*(u-mu)/s2)-1;
      let a=A[i]+dt*g; a=a<0?0:a>1?1:a; A[i]=a; mass+=a; }
    return mass;
  }
  return {setKernel, step, NC};
}


let N,L,A,U,mu,sigma,dt,gen=0;
function send(ms){
  let mass=0; for(let i=0;i<A.length;i++) mass+=A[i];
  const hist=new Float32Array(200);
  for(let i=0;i<U.length;i++){ const u=U[i]; if(u>0.002){ const b=Math.min(199,(u*200)|0); hist[b]++; } }
  const copy=A.slice();
  postMessage({type:'state',A:copy,mass,ms,gen,hist},[copy.buffer,hist.buffer]);
}
onmessage=e=>{ const d=e.data;
  if(d.type==='init'){ N=d.N; L=makeLenia(N); A=new Float32Array(d.A); U=new Float32Array(N*N*N); L.setKernel(d.R,d.beta); mu=d.mu; sigma=d.sigma; dt=d.dt; gen=0; send(0); }
  else if(d.type==='kernel'){ L.setKernel(d.R,d.beta); }
  else if(d.type==='growth'){ mu=d.mu; sigma=d.sigma; dt=d.dt; }
  else if(d.type==='state'){ A=new Float32Array(d.A); if(d.resetGen) gen=0; send(0); }
  else if(d.type==='step'){ const t0=performance.now(); for(let k=0;k<d.k;k++) L.step(A,U,mu,sigma,dt); gen+=d.k; send((performance.now()-t0)/d.k); }
};
