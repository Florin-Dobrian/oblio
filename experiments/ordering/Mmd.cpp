// Mmd.cc — Multiple Minimum Degree ordering.
//
// Ported from oblio 0.9, which ported from Liu/Sparspak genmmd (1985).
// All functions in anonymous namespace; 1-based indexing convention throughout.
// Single entry point: mmd_order (0-based interface, matches AMD convention).

#include <cstdint>
#include <vector>

namespace {

// ============================================================================
// genmmd and helpers — 1-based indexing, ported from Sparspak via oblio 0.9.
// ============================================================================

static void mmdelm(int,int[],int[],int[],int[],int[],int[],int[],int[],int,int);
static int  mmdint(int,int[],int[],int[],int[],int[],int[],int[]);
static void mmdnum(int,int[],int[],int[]);
static void mmdupd(int,int,int[],int[],int,int*,int[],int[],int[],int[],int[],int[],int,int*);

static void genmmd(int neqns,int xadj[],int adjncy[],int invp[],int perm[],
                   int delta,int head[],int qsize[],int list[],int marker[],
                   int maxint,int*ncsub)
{
    if(neqns<=0)return;
    xadj--;adjncy--;invp--;perm--;head--;qsize--;list--;marker--;
    *ncsub=0;
    mmdint(neqns,xadj,head,invp,perm,qsize,list,marker);
    int num=1;
    int nextmd=head[1];
    while(nextmd>0){int mn=nextmd;nextmd=invp[mn];marker[mn]=maxint;invp[mn]=-num;num++;}
    if(num>neqns)goto n1000;
    {int tag=1,mdeg=2;head[1]=0;
    while(1){
        while(head[mdeg]<=0)mdeg++;
        int mdlmt=mdeg+delta,ehead=0;
n500:   {int mn=head[mdeg];
        while(mn<=0){mdeg++;if(mdeg>mdlmt)goto n900;mn=head[mdeg];}
        int nx=invp[mn];head[mdeg]=nx;if(nx>0)perm[nx]=-mdeg;
        invp[mn]=-num;*ncsub+=mdeg+qsize[mn]-2;
        if((num+qsize[mn])>neqns)goto n1000;
        tag++;if(tag>=maxint){tag=1;for(int i=1;i<=neqns;i++)if(marker[i]<maxint)marker[i]=0;}
        mmdelm(mn,xadj,adjncy,head,invp,perm,qsize,list,marker,maxint,tag);
        num+=qsize[mn];list[mn]=ehead;ehead=mn;
        if(delta>=0)goto n500;}
n900:   if(num>neqns)goto n1000;
        mmdupd(ehead,neqns,xadj,adjncy,delta,&mdeg,head,invp,perm,qsize,list,marker,maxint,&tag);
    }}
n1000:mmdnum(neqns,perm,invp,qsize);
    xadj++;adjncy++;invp++;perm++;head++;qsize++;list++;marker++;
}

static void mmdelm(int md,int xadj[],int adjncy[],int head[],
                   int fwd[],int bwd[],int qsize[],int list[],
                   int marker[],int maxint,int tag)
{
    marker[md]=tag;
    int is=xadj[md],it=xadj[md+1]-1,el=0,rl=is,rm=it;
    for(int i=is;i<=it;i++){
        int nb=adjncy[i];if(nb==0)break;
        if(marker[nb]<tag){marker[nb]=tag;
            if(fwd[nb]<0){list[nb]=el;el=nb;}
            else{adjncy[rl]=nb;rl++;}}}
    while(el>0){
        adjncy[rm]=-el;int lk=el;
n400:   int js=xadj[lk],jt=xadj[lk+1]-1;
        for(int j=js;j<=jt;j++){
            int nd=adjncy[j];lk=-nd;
            if(nd<0)goto n400;if(nd==0)break;
            if(marker[nd]<tag&&fwd[nd]>=0){
                marker[nd]=tag;
                while(rl>=rm){lk=-adjncy[rm];rl=xadj[lk];rm=xadj[lk+1]-1;}
                adjncy[rl]=nd;rl++;}}
        el=list[el];}
    if(rl<=rm)adjncy[rl]=0;
    int lk=md;
n1100:{int js=xadj[lk],jt=xadj[lk+1]-1;
    for(int i=js;i<=jt;i++){
        int rn=adjncy[i];lk=-rn;
        if(rn<0)goto n1100;if(rn==0)return;
        int pv=bwd[rn];
        if(pv!=0&&pv!=(-maxint)){
            int nx=fwd[rn];if(nx>0)bwd[nx]=pv;if(pv>0)fwd[pv]=nx;
            if(pv<0)head[-pv]=nx;}
        int js2=xadj[rn],jt2=xadj[rn+1]-1,xq=js2;
        for(int j=js2;j<=jt2;j++){int nb=adjncy[j];if(nb==0)break;if(marker[nb]<tag){adjncy[xq]=nb;xq++;}}
        int nq=xq-js2;
        if(nq<=0){qsize[md]+=qsize[rn];qsize[rn]=0;marker[rn]=maxint;fwd[rn]=-md;bwd[rn]=-maxint;}
        else{fwd[rn]=nq+1;bwd[rn]=0;adjncy[xq]=md;xq++;if(xq<=jt2)adjncy[xq]=0;}
    }}
}

static int mmdint(int neqns,int xadj[],int head[],int fwd[],int bwd[],int qsize[],int list[],int marker[])
{
    for(int nd=1;nd<=neqns;nd++){head[nd]=0;qsize[nd]=1;marker[nd]=0;list[nd]=0;}
    for(int nd=1;nd<=neqns;nd++){
        int dg=xadj[nd+1]-xadj[nd];if(dg==0)dg=1;
        int fn=head[dg];fwd[nd]=fn;head[dg]=nd;if(fn>0)bwd[fn]=nd;bwd[nd]=-dg;}
    return 0;
}

static void mmdnum(int neqns,int perm[],int invp[],int qsize[])
{
    for(int nd=1;nd<=neqns;nd++){if(qsize[nd]<=0)perm[nd]=invp[nd];else perm[nd]=-invp[nd];}
    for(int nd=1;nd<=neqns;nd++){
        if(perm[nd]<=0){
            int fa=nd;while(perm[fa]<=0)fa=-perm[fa];
            int rt=fa,nm=perm[rt]+1;invp[nd]=-nm;perm[rt]=nm;
            fa=nd;int nf=-perm[fa];while(nf>0){perm[fa]=-rt;fa=nf;nf=-perm[fa];}}}
    for(int nd=1;nd<=neqns;nd++){int nm=-invp[nd];invp[nd]=nm;perm[nm]=nd;}
}

static void mmdupd(int ehead,int neqns,int xadj[],int adjncy[],int delta,int*mdeg,
                   int head[],int fwd[],int bwd[],int qsize[],int list[],int marker[],int maxint,int*tag)
{
    /* All variables hoisted to avoid goto-crosses-initialization errors. */
    int md0=*mdeg+delta,el=ehead;
    int mt,q2h,qxh,dg0,lk,en,iq2,dg;
    int is,it,nb,is2,it2,nd,js,jt,fn;
    (void)fn;
n100:
    if(el<=0)return;
    mt=*tag+md0;
    if(mt>=maxint){*tag=1;for(int i=1;i<=neqns;i++)if(marker[i]<maxint)marker[i]=0;mt=*tag+md0;}
    q2h=0;qxh=0;dg0=0;lk=el;
n400:
    is=xadj[lk];it=xadj[lk+1]-1;
    for(int i=is;i<=it;i++){
        nb=adjncy[i];lk=-nb;if(nb<0)goto n400;if(nb==0)break;
        if(qsize[nb]!=0){dg0+=qsize[nb];marker[nb]=mt;
            if(bwd[nb]==0){if(fwd[nb]!=2){list[nb]=qxh;qxh=nb;}else{list[nb]=q2h;q2h=nb;}}}}
    en=q2h;iq2=1;
n900:
    if(en<=0)goto n1500;if(bwd[en]!=0)goto n2200;
    (*tag)++;dg=dg0;
    is=xadj[en];nb=adjncy[is];if(nb==el)nb=adjncy[is+1];lk=nb;
    if(fwd[nb]>=0){dg+=qsize[nb];goto n2100;}
n1000:
    is2=xadj[lk];it2=xadj[lk+1]-1;
    for(int i=is2;i<=it2;i++){
        nd=adjncy[i];lk=-nd;
        if(nd!=en){if(nd<0)goto n1000;if(nd==0)goto n2100;
        if(qsize[nd]!=0){if(marker[nd]<*tag){marker[nd]=*tag;dg+=qsize[nd];}
        else if(bwd[nd]==0){if(fwd[nd]==2){qsize[en]+=qsize[nd];qsize[nd]=0;marker[nd]=maxint;fwd[nd]=-en;bwd[nd]=-maxint;}
        else if(bwd[nd]==0)bwd[nd]=-maxint;}}}}
    goto n2100;
n1500:
    en=qxh;iq2=0;
n1600:
    if(en<=0)goto n2300;if(bwd[en]!=0)goto n2200;
    (*tag)++;dg=dg0;
    is=xadj[en];it=xadj[en+1]-1;
    for(int i=is;i<=it;i++){nb=adjncy[i];if(nb==0)break;
        if(marker[nb]<*tag){marker[nb]=*tag;lk=nb;
            if(fwd[nb]>=0){dg+=qsize[nb];}
            else{
n1700:
                js=xadj[lk];jt=xadj[lk+1]-1;
                for(int j=js;j<=jt;j++){nd=adjncy[j];lk=-nd;if(nd<0)goto n1700;if(nd==0)break;
                    if(marker[nd]<*tag){marker[nd]=*tag;dg+=qsize[nd];}}}}}
n2100:
    dg=dg-qsize[en]+1;if(dg<1)dg=1;
    fn=head[dg];
    fwd[en]=fn;bwd[en]=-dg;if(fn>0)bwd[fn]=en;head[dg]=en;if(dg<*mdeg)*mdeg=dg;
n2200:
    en=list[en];if(iq2==1)goto n900;goto n1600;
n2300:
    *tag=mt;el=list[el];goto n100;
}

} // anonymous namespace

// ============================================================================
// Public entry point — 0-based interface.
//
// Input:  n, colPtr[n+1], rowIdx[nnz] — off-diagonal-only symmetric CSC.
// Output: perm[n] — new-to-old permutation (P[k] = original node of pivot k).
//         invp[n] — old-to-new permutation.
// ============================================================================

void mmd_order(int n, const int colPtr[], const int rowIdx[],
               int perm[], int invp[])
{
    if (n <= 0) return;
    int nnz = colPtr[n];

    // Convert 0-based to 1-based for genmmd.
    std::vector<int> xadj(n + 2), adj(nnz + 1);
    std::vector<int> head(n + 1), qsize(n + 1), list(n + 1), marker(n + 1);

    for (int j = 0; j <= n; j++) xadj[j + 1] = colPtr[j] + 1;
    for (int k = 0; k < nnz; k++) adj[k + 1] = rowIdx[k] + 1;

    // invp and perm passed as 1-based (offset by 1)
    std::vector<int> invp1(n + 1), perm1(n + 1);

    int ncsub = 0;
    genmmd(n, xadj.data() + 1, adj.data() + 1,
           invp1.data() + 1, perm1.data() + 1,
           0, head.data() + 1, qsize.data() + 1,
           list.data() + 1, marker.data() + 1,
           8388607, &ncsub);

    // Convert back to 0-based.
    for (int j = 0; j < n; j++) {
        invp[j] = invp1[j + 1] - 1;  // old-to-new
        perm[j] = perm1[j + 1] - 1;  // new-to-old
    }
}
