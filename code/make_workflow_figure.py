"""Draw Figure 1 using Matplotlib primitives, with no image-model or image input.

The default canvas is 7.4 inches wide. The smallest body label is 8 points
(about 7.4 points when printed at 6.8 inches), larger than the 6-point caption
in the supplied OUP proof. Text extents are checked against their panels.
"""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle


def build(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.4, 4.7), facecolor='white')
    ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, 100), ylim=(0, 64)); ax.axis('off')
    ink = '#17243F'; navy = '#173D70'
    panels = [(1, 2, 23, 53), (26, 2, 23, 53), (51, 2, 23, 53), (76, 2, 23, 53)]
    colors = ['#65717E', '#2E78C4', '#7042C1', '#2E914F']
    fills = ['#F3F6FA', '#EAF4FC', '#F5EEFD', '#ECF8EF']
    checks = []
    def box(b, fill, edge):
        x,y,w,h=b
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.15,rounding_size=1.2',facecolor=fill,edgecolor=edge,lw=.85))
    def text(x,y,s,size=9,weight='normal',color=ink,ha='center',panel=None,va='center'):
        t=ax.text(x,y,s,fontsize=size,fontweight=weight,color=color,ha=ha,va=va,linespacing=1.14)
        if panel is not None: checks.append((t,panels[panel]))
        return t
    text(50,60,'GoCo workflow',16,'bold',navy)
    titles=['Predictor\noutputs','Candidate\npolicies','Certification\non C','Common\nevaluation on E']
    for j,(b,c,f,title) in enumerate(zip(panels,colors,fills,titles)):
        box(b,f,c);x,y,w,h=b
        ax.add_patch(Circle((x+2.8,51),1.5,color=c))
        text(x+2.8,51,str(j+1),9,'bold','white')
        text(x+13,51,title,10,'bold',panel=j)
    text(12.5,44,'gene-GO scores',10,panel=0)
    text(12.5,40.8,r'$s_{ig}$',10,panel=0)
    text(12.5,36,'sources behind each\ncandidate call',9,panel=0)
    ramp=['#DCEAF5','#BDD8EC','#8BB9DB','#5799CA','#2E78C4']
    for r in range(4):
        for c in range(5):
            ax.add_patch(Rectangle((7.6+c*2,23+r*2),2,2,facecolor=ramp[(r+c)%5],edgecolor='white',lw=.4))
    text(12.5,21,'GO terms',8,panel=0)
    t=text(5.4,27,'genes',8,panel=0);t.set_rotation(90)
    box((2.8,4,19.4,13.2),'#E9EEF3','#BCC8D4')
    text(12.5,15.2,'Data split',10,'bold',navy,panel=0)
    text(12.5,11,'ranking T: 10%\ncertification C: 60%\nevaluation E: 30%',8.5,panel=0)
    text(12.5,5.7,'uniform random\npartition',8.2,'bold',panel=0)
    text(37.5,43.3,'Nested paths using\npredictor-native\nstructure',8.0,'bold',navy,panel=1)
    box((27.6,26,19.8,12.7),'#FDEEEF','#F06B75')
    text(37.5,35.8,'GoCo-M',11,'bold','#D23B4B',panel=1)
    text(37.5,30,'shared module evidence\npartial admission of\ngenes within grid steps',8.3,panel=1)
    box((27.6,10.8,19.8,13),'#FFF3E6','#F28A1D')
    text(37.5,20.9,'GoCo-N',11,'bold','#D96A00',panel=1)
    text(37.5,15.3,'gene-specific evidence\nfrom neighbourhoods\nindividual-call ordering',8.0,panel=1)
    text(37.5,6.2,'Paths fixed using T\nbefore certification',8.5,panel=1)
    text(62.5,43,'Test candidate policies\nin path order',9,panel=2)
    text(62.5,36,r'$\pi_1,\ \pi_2,\ \ldots,\ \pi_K$',11,color=navy,panel=2)
    text(62.5,29,'Fixed-sequence\nLearn-Then-Test',10,'bold',navy,panel=2)
    text(62.5,22,r'Stop at first $p_k>\delta$',9,panel=2)
    box((52.5,5,20,11.7),'#E7DAFA','#E7DAFA')
    text(62.5,14.5,'Proposition 1',9.5,'bold',navy,panel=2)
    text(62.5,10.4,r'$\Pr\{R(\widehat\pi)>\alpha\}\leq\delta$'+'\nasymptotically',8.5,color=navy,panel=2)
    text(62.5,5.9,'under stated assumptions',8,color=navy,panel=2)
    text(87.5,43,'Held-out gene-level\nTruePath FDP',9.5,panel=3)
    text(87.5,33,'Reference-supported\nreleased calls',9.5,panel=3)
    # Schematic metric symbols, not numerical results.
    for k,h in enumerate([6,4.5,3,2]):
        ax.add_patch(Rectangle((81.7+k*3.3,18),2.2,h,facecolor='#4CA36A',edgecolor='none'))
    text(87.5,15.8,'yield and risk diagnostics',8,panel=3)
    text(87.5,8,'Same frozen reference\nfor every method',9,'bold',panel=3)
    for x in (24.2,49.2,74.2):
        ax.add_patch(FancyArrowPatch((x,31),(x+1.5,31),arrowstyle='-|>',mutation_scale=13,color='#65717E',lw=.9))
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for t,(x,y,w,h) in checks:
        b=t.get_window_extent(renderer).transformed(ax.transData.inverted())
        if b.x0<x+.25 or b.x1>x+w-.25 or b.y0<y+.25 or b.y1>y+h-.25:
            raise RuntimeError('Text crosses panel: '+t.get_text())
    fig.savefig(outdir/'GoCo_Workflow.png',dpi=300,facecolor='white')
    fig.savefig(outdir/'GoCo_Workflow.pdf',facecolor='white')
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=Path(__file__).resolve().parents[1]/'paper'/'figures')
    build(p.parse_args().outdir)
