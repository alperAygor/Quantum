import os  # CIKTI_DIR: figurlerin yazilacagi klasor (varsayilan: calisilan klasor)
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
BLUE='#2a78d6'; ORANGE='#eb6834'; INK='#0b0b0b'; INK2='#52514e'; GRID='#e4e3df'; SURF='#fcfcfb'
d=np.loadtxt('profil.xvg',comments=['#','@']); x,w=d[:,0],d[:,1]
m=(x>-1.92)&(x<1.92); x,w=x[m],w[m]
lo=w[(x>=-1.9)&(x<=-1.0)].mean(); hi=w[(x>=1.0)&(x<=1.9)].mean(); ref=(lo+hi)/2
W=w-ref
h=np.loadtxt('histo.xvg',comments=['#','@'])
plt.rcParams.update({'font.size':11,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2})
fig,(a1,a2)=plt.subplots(2,1,figsize=(10,8.2),gridspec_kw={'height_ratios':[3,1.3]},sharex=True,facecolor=SURF)
for a in (a1,a2):
    a.set_facecolor(SURF); a.grid(axis='y',color=GRID,lw=0.8); a.spines[['top','right']].set_visible(False)
a1.plot(x,W,color=BLUE,lw=2,solid_capstyle='round')
a1.axhline(lo-ref,color=INK2,lw=1,ls=(0,(4,3))); a1.axhline(hi-ref,color=INK2,lw=1,ls=(0,(4,3)))
a1.text(-1.88,lo-ref-1.5,f'Su (dar kenar tarafı): {lo-ref:+.1f} kJ/mol',color=INK2,fontsize=10)
a1.text(1.88,hi-ref+0.6,f'Su (geniş kenar tarafı): {hi-ref:+.1f} kJ/mol',color=INK2,fontsize=10,ha='right')
i0=np.argmin(W); a1.plot(x[i0],W[i0],'o',ms=8,color=BLUE,mec=SURF,mew=2)
a1.annotate(f'Bağlı durum (kavite merkezi)\nW = {W[i0]:.1f} kJ/mol, ξ = {x[i0]:+.2f} nm',xy=(x[i0],W[i0]),xytext=(0.35,-17.6),color=INK,fontsize=10,arrowprops=dict(arrowstyle='-',color=INK2,lw=1))
j=np.argmin(abs(x+0.6))
a1.annotate('Dar kenarda ara durum\n(ξ ≈ −0.45 … −0.75 nm)',xy=(x[j],W[j]),xytext=(-1.75,-9.5),color=INK,fontsize=10,arrowprops=dict(arrowstyle='-',color=INK2,lw=1))
k=np.argmin(abs(x-0.05))
a1.annotate('ξ≈0\'da basamak:\niki çekmenin birleşim noktası',xy=(x[k],W[k]),xytext=(0.35,-13.2),color=INK,fontsize=10,arrowprops=dict(arrowstyle='-',color=INK2,lw=1))
a1.axvline(0,color=ORANGE,lw=1,ls=(0,(2,2)))
a1.text(0.03,6.6,'ξ<0 pencereleri: dar kenar çekmesinden  |  ξ≥0 pencereleri: geniş kenar çekmesinden',color=INK2,fontsize=9,va='top')
a1.set_ylabel('W(ξ)  [kJ/mol]'); a1.set_ylim(-19,7)
a1.set_title('β-CD + siklopentanol: ortalama kuvvet potansiyeli (PMF)',loc='left',color=INK,fontsize=13,pad=18)
a1.text(0,1.01,'Umbrella sampling, 37 pencere × 2 ns, WHAM, 300 K. Sıfır: iki su platosunun ortalaması.',transform=a1.transAxes,color=INK2,fontsize=9.5,va='bottom')
for c in range(1,h.shape[1]):
    a2.plot(h[:,0],h[:,c],color=BLUE,lw=1,alpha=0.75)
a2.set_ylabel('Örnek sayısı'); a2.set_xlabel('ξ: misafirin CD ekseni boyunca konumu [nm]   (− dar kenar,  + geniş kenar)')
a2.set_title('Pencere histogramları (örtüşme kontrolü)',loc='left',color=INK,fontsize=11)
a2.set_xlim(-1.95,1.95)
fig.tight_layout(); fig.savefig(os.environ.get("CIKTI_DIR", ".") + '/PMF_bCD_siklopentanol.png',dpi=160,facecolor=SURF)
