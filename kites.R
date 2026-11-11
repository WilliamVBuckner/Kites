library(ape); library(phytools); library(phangorn); library(brms); library(maps); library(ggplot2); library(magick)
set.seed(1)
options(mc.cores=4)

dir.create("cache",showWarnings=FALSE)
key <- function(...) {f <- tempfile(); saveRDS(list(...),f); substr(unname(tools::md5sum(f)),1,8)}
get <- function(url,f) {f <- file.path("cache",paste0(key(url),"_",f)); if(!file.exists(f)) {download.file(url,paste0(f,".part"),mode="wb"); file.rename(paste0(f,".part"),f)}; f}
gray <- "https://raw.githubusercontent.com/D-PLACE/dplace-data/v2.2.1/phylogenies/gray_et_al2009/"
eau <- "https://raw.githubusercontent.com/D-PLACE/dplace-dataset-ea/v3.3/cldf/"
gl <- read.csv(get("https://raw.githubusercontent.com/glottolog/glottolog-cldf/v5.3/cldf/values.csv","values.csv"))
gl <- gl[gl$Parameter_ID=="classification",]
lin <- setNames(strsplit(paste(gl$Value,gl$Language_ID,sep="/"),"/"),gl$Language_ID); rm(gl)
tx <- read.csv(get(paste0(gray,"taxa.csv"),"taxa.csv")); tx <- tx[tx$glottocode!="",]
pt <- read.nexus(get(paste0(gray,"posterior.trees"),"posterior.trees")); hp <- key(pt)
st <- read.nexus(get(paste0(gray,"summary.trees"),"summary.trees"))
ea <- read.csv(get(paste0(eau,"data.csv"),"ea.csv"))
eas <- read.csv(get(paste0(eau,"societies.csv"),"eas.csv"))
eav <- function(v,id) {x <- ea[ea$Var_ID==v,]; as.numeric(sub(".*-","",x$Code_ID[match(id,x$Soc_ID)]))}
gz <- get("https://github.com/rbouckaert/global-language-tree-pipeline/releases/download/v2.0.0/global-language-trees-6636-taxa.tgz","global.tgz")
gf <- file.path(sub("\\.tgz$","",gz),"global-language-trees-6636-taxa.trees")
if(!file.exists(gf)) untar(gz,exdir=dirname(gf))

# society level
d0 <- read.csv("kites.csv")
table(d0$source[d0$relevant==1]); length(unique(d0$society[d0$relevant==1]))
d <- d0[d0$relevant==1 & d0$glottocode!="",]
dom <- c("play","contest","gambling","fishing","ritual","divination","festival","myth","sound","utility","prohibition","commerce")
fl <- strsplit(d$fliers,"; ")
d$adult <- sapply(fl,function(x) any(x %in% c("men","women","adults","elders","chiefs","priests","everyone")))
d$child <- sapply(fl,function(x) any(x %in% c("boys","girls","children","youths")))
d$deny_kite <- d$denial==1 & grepl("^kites$|^kite-flying|play kites",d$denied)
d$deny_fish <- d$denial==1 & grepl("fishing",d$denied) & !grepl("disuse",d$denied)
d$regional <- d$scope=="regional"
d$aff <- rowSums(d[,c(dom,"construction","other")])>0 & !d$deny_kite
d$made <- d$aff & d$narrative==0 & rowSums(d[,c(setdiff(dom,"myth"),"construction")])>0
d$ritual2 <- ifelse(d$pid %in% c("cha_15","cha_16","cha_18"),0,d$ritual)
first <- function(x) {x <- x[!is.na(x) & x!=""]; if(length(x)) x[1] else NA}
s <- do.call(rbind,lapply(split(d,d$society),function(x) data.frame(
  society=x$society[1],glottocode=x$glottocode[1],lineage=x$lineage[1],ea_id=first(x$ea_id),
  lat=as.numeric(first(x$latitude)),lon=as.numeric(first(x$longitude)),
  ehraf=any(x$source=="eHRAF"),lit=any(x$source!="eHRAF"),present=any(x$aff),made=any(x$made),
  t(sapply(dom,function(v) max(x[[v]]))),ritual2=max(x$ritual2),
  deny_kite=any(x$deny_kite),deny_fish=any(x$deny_fish),deny_fish_l=any(x$deny_fish & !x$regional),deny_kite_l=any(x$deny_kite & !x$regional),adult=any(x$adult),child=any(x$child))))
s$family <- sapply(strsplit(s$lineage,"/"),`[`,1)
s$fd <- eav("EA003",s$ea_id); s$jh <- eav("EA033",s$ea_id)
p <- s[s$present,]
nrow(s); nrow(p); table(p$family=="aust1307"); table(ehraf=p$ehraf,lit=p$lit)
cbind(n=colSums(p[,dom]),prop=round(colMeans(p[,dom]),2))
table(adult=p$adult,child=p$child)

# frames, grafting and reconstruction
a <- s[s$family=="aust1307",]
a$poly <- grepl("poly1242",a$lineage)
fk <- function(a,regional=TRUE) {df <- if(regional) a$deny_fish else a$deny_fish_l; dk <- if(regional) a$deny_kite else a$deny_kite_l
  ifelse(a$fishing==1 & (df|dk),0.5,ifelse(a$fishing==1,1,ifelse(df|dk|(regional & a$poly),0,NA)))}
kk <- function(a,v=a$made) ifelse(v & a$deny_kite,0.5,ifelse(v,1,ifelse(a$deny_kite,0,NA)))
cp <- function(x,y) {n <- min(length(x),length(y)); m <- which(x[1:n]!=y[1:n]); if(length(m)) m[1]-1 else n}
match_tips <- function(u,tx) {
  u$taxon <- tx$taxon[match(u$glottocode,tx$glottocode)]
  u$graft <- sapply(seq_len(nrow(u)),function(i) {
    if(!is.na(u$taxon[i])) return(u$taxon[i])
    k <- sapply(lin[tx$glottocode],cp,x=strsplit(u$lineage[i],"/")[[1]])
    if(max(k)==0) "" else paste(tx$taxon[k==max(k)],collapse=";")})
  u
}
code <- function(a,f=fk(a),tt=tx,conf=0.5) {
  a$fish <- f; a <- a[!is.na(f),]
  u <- aggregate(cbind(fish,lat,lon,fd,ritual,ritual2,play)~glottocode+lineage,a,function(x) if(all(x %in% 0:1) & length(unique(x))>1) 0.5 else mean(x,na.rm=TRUE),na.action=na.pass)
  u$fish[!u$fish %in% c(0,1)] <- conf; u$ritual <- ceiling(u$ritual); u$ritual2 <- ceiling(u$ritual2); u$play <- ceiling(u$play)
  match_tips(u,tt)
}
graft <- function(t,u) {
  t <- keep.tip(t,intersect(unique(c(u$taxon,unlist(strsplit(u$graft,";")))),t$tip.label))
  for(i in which(is.na(u$taxon))) {
    g <- strsplit(u$graft[i],";")[[1]]; g <- g[g %in% t$tip.label]
    h <- node.depth.edgelength(t); top <- max(h)
    if(length(g)==1) {w <- which(t$tip.label==g); e <- t$edge.length[t$edge[,2]==w]/2; t <- bind.tip(t,u$glottocode[i],edge.length=top-h[w]+e,where=w,position=e)}
    else {w <- if(length(g)) getMRCA(t,g) else Ntip(t)+1; t <- bind.tip(t,u$glottocode[i],edge.length=top-h[w],where=w,position=0)}
  }
  k <- !is.na(u$taxon) & u$taxon %in% t$tip.label
  t$tip.label[match(u$taxon[k],t$tip.label)] <- u$glottocode[k]
  multi2di(keep.tip(t,u$glottocode))
}
cl <- c(malayo_polynesian="mala1545",oceanic="ocea1241",micronesian="micr1243",central_pacific="cent2060",polynesian="poly1242")
asr <- function(u,n=100,trees=pt,pi="fitzjohn") {
  x <- cbind("0"=1-u$fish,"1"=u$fish); rownames(x) <- u$glottocode; set.seed(1)
  r <- do.call(rbind,lapply(sample(length(trees),n),function(i) {
    t <- graft(trees[[i]],u)
    f <- lapply(c("ER","ARD"),function(m) fitMk(t,x,model=m,pi=pi))
    w <- exp(-0.5*(sapply(f,AIC)-min(sapply(f,AIC)))); w <- w/sum(w)
    an <- sapply(f,function(z) ancr(z)$ace[,2]); q <- f[[2]]$rates[f[[2]]$index.matrix]; dim(q) <- c(2,2)
    nd <- sapply(cl,function(k) getMRCA(t,u$glottocode[grepl(k,u$lineage)]))
    h <- max(node.depth.edgelength(t))-node.depth.edgelength(t)
    data.frame(tree=i,node=names(nd),p=as.vector(an[as.character(nd),] %*% w),w_er=w[1],q_er=f[[1]]$rates,gain=q[1,2],loss=q[2,1],age=h[nd])
  }))
  r$node <- factor(r$node,names(cl)); r
}
asr3 <- function(u,n=100,trees=pt) {
  x <- as.matrix(u[,c("s0","s1","s2")])+0; dimnames(x) <- list(u$glottocode,0:2); set.seed(1)
  m <- lapply(list(c(0,1,0,1,0,1,0,1,0),c(0,1,0,2,0,3,0,4,0)),matrix,3,3,dimnames=list(0:2,0:2))
  r <- do.call(rbind,lapply(sample(length(trees),n),function(i) {
    t <- graft(trees[[i]],u)
    f <- lapply(m,function(k) fitMk(t,x,model=k,pi="fitzjohn"))
    w <- exp(-0.5*(sapply(f,AIC)-min(sapply(f,AIC)))); w <- w/sum(w)
    an <- lapply(f,function(z) ancr(z)$ace)
    nd <- as.character(sapply(cl,function(k) getMRCA(t,u$glottocode[grepl(k,u$lineage)])))
    h <- max(node.depth.edgelength(t))-node.depth.edgelength(t)
    data.frame(tree=i,node=names(cl),p=as.vector(1-sapply(an,function(z) z[nd,1]) %*% w),pf=as.vector(sapply(an,function(z) z[nd,3]) %*% w),w_er=w[1],q_er=NA,gain=NA,loss=NA,age=h[as.integer(nd)])
  }))
  r$node <- factor(r$node,names(cl)); r
}
joint <- function(k,f) {j <- merge(k[,c("glottocode","fish")],f[,c("glottocode","fish")],by="glottocode",all=TRUE)
  j <- cbind(j,rbind(k,f)[match(j$glottocode,c(k$glottocode,f$glottocode)),c("lineage","taxon","graft")])
  a <- j$fish.x; b <- j$fish.y
  j$s0 <- (is.na(a)|a<1)&(is.na(b)|b<1); j$s1 <- (is.na(a)|a>0)&(is.na(b)|b<1); j$s2 <- (is.na(a)|a>0)&(is.na(b)|b>0)
  j$fish <- j$s0+2*j$s1+4*j$s2; j}
sm <- function(r) aggregate(cbind(p,age)~node,r,function(z) round(quantile(z,c(.5,.1,.9)),4))
cached <- function(f,k,x) {f <- file.path("cache",paste0(f,"_",k,".rds")); if(file.exists(f)) readRDS(f) else {saveRDS(x,f); x}}
fit <- function(f,x,...) cached(f,key(...,deparse(substitute(x))),x)
res <- list()
run <- function(f,u,n=100,trees=pt,tr=hp,fn=asr,...) {r <- cached(f,key(u[,c("glottocode","lineage","fish","graft")],n,tr,deparse(fn),deparse(graft),cl,...),fn(u,n,trees,...)); res[[f]] <<- r; r}
smf <- function(r) aggregate(pf~node,r,function(z) round(quantile(z,c(.5,.1,.9)),4))

# fishing kite on the austronesian trees
u <- code(a)
table(u$fish); nrow(u); sum(is.na(u$taxon))
r0 <- run("asr",u); sm(r0)
summary(r0$w_er[!duplicated(r0$tree)])
u1 <- code(a,fk(a,regional=FALSE)); table(u1$fish); sm(run("asr_noregional",u1,50))
u2 <- code(a,conf=1); table(u2$fish); sm(run("asr_conflicts",u2,50))
u3 <- code(a,conf=0); table(u3$fish); sm(run("asr_conflicts_absent",u3,50))
sm(run("asr_oceanic",u[grepl("ocea1241",u$lineage),],50))
sm(run("asr_equalroot",u,50,pi="equal"))
sm(run("asr_matched",u[!is.na(u$taxon),],50))

# kites themselves: made or flown, then any record
uk <- code(a,kk(a))
table(uk$fish); nrow(uk); sum(is.na(uk$taxon))
rk <- run("asr_kite",uk); sm(rk)
sm(run("asr_kite_absent",code(a,kk(a),conf=0),50))
sm(run("asr_kite_present",code(a,kk(a),conf=1),50))
ud <- code(a,kk(a,a$present)); table(ud$fish)
rd <- run("asr_kite_any",ud); sm(rd)
uj <- joint(uk,u); table(uj$fish)
rj <- run("asr_joint",uj,fn=asr3); sm(rj); smf(rj); res$asr_joint_fishing <- transform(rj,p=pf)

# the same on the global trees of bouckaert et al.
rm(pt); invisible(gc())
set.seed(2)
gl2 <- readLines(gf,warn=FALSE)
i0 <- grep("Taxlabels",gl2); glab <- trimws(gl2[(i0+1):(grep("^\\s*;",gl2)[grep("^\\s*;",gl2)>i0][1]-1)])
gpt <- lapply(sample(grep("^tree ",gl2),100),function(i) {t <- read.tree(text=sub("^tree [^=]+= *","",gl2[i])); t$tip.label <- glab[as.integer(t$tip.label)]; t})
rm(gl2); invisible(gc())
gx <- data.frame(taxon=glab,glottocode=glab); hg <- key(gpt)
ug <- code(a,tt=gx); sum(is.na(ug$taxon))
rg <- run("asr_global",ug,100,gpt,hg); sm(rg)
rgk <- run("asr_global_kite",code(a,kk(a),tt=gx),100,gpt,hg); sm(rgk)
sm(run("asr_global_kite_any",code(a,kk(a,a$present),tt=gx),100,gpt,hg))
rgj <- run("asr_global_joint",joint(code(a,kk(a),tt=gx),ug),100,gpt,hg,asr3); sm(rgj); smf(rgj); res$asr_global_joint_fishing <- transform(rgj,p=pf)

# descent or proximity
xyz <- function(lat,lon) {la <- lat*pi/180; lo <- lon*pi/180; cbind(c1=cos(la)*cos(lo),c2=cos(la)*sin(lo),c3=sin(la))}
set.seed(1)
b <- u[u$fish %in% 0:1,]; tb <- graft(st,b)
A <- vcv.phylo(tb,corr=TRUE)
b <- cbind(b,xyz(b$lat,b$lon))
share <- function(m) {v <- as_draws_df(m); ph <- v[[grep("^sd_.*Intercept",names(v))]]^2; sp <- v[[grep("^sdgp",names(v))]]^2
  round(quantile(ph/(ph+sp),c(.5,.05,.95)),2)}
m_fish <- fit("m_fish",brm(fish~1+(1|gr(glottocode,cov=A))+gp(c1,c2,c3),data=b,data2=list(A=A),family=bernoulli,seed=1,control=list(adapt_delta=.99)),b[,c("fish","glottocode","c1","c2","c3")],A)
summary(m_fish); share(m_fish)

# fishing dependence and the sacred kite
eff <- function(m,v) {z <- as_draws_df(m)[[paste0("b_",v)]]; round(c(est=mean(z),quantile(z,c(.025,.975)),pd=max(mean(z>0),mean(z<0))),2)}
b$y <- b$fish; b$id <- b$glottocode; b$x <- b$ritual
bf <- b[!is.na(b$fd),]; bf$x <- as.vector(scale(bf$fd))
Af <- vcv.phylo(keep.tip(tb,bf$id),corr=TRUE)
m_fd <- fit("m_fd",brm(y~x+(1|gr(id,cov=V))+gp(c1,c2,c3),data=bf,data2=list(V=Af),family=bernoulli,prior=prior(normal(0,2),class=b),seed=1,control=list(adapt_delta=.99)),bf[,c("y","x","id","c1","c2","c3")],Af)
m_rit <- fit("m_rit",update(m_fd,newdata=b,data2=list(V=A),seed=1,control=list(adapt_delta=.99)),b[,c("y","x","id","c1","c2","c3")],A,stancode(m_fd))
nrow(bf); table(fish=bf$fish,fd=bf$fd); eff(m_fd,"x")
table(fish=b$fish,ritual=b$ritual); table(fish=b$fish,play=b$play); eff(m_rit,"x")
b$x <- b$ritual2; m_rit2 <- fit("m_rit2",update(m_fd,newdata=b,data2=list(V=A),seed=1,control=list(adapt_delta=.99)),b[,c("y","x","id","c1","c2","c3")],A,stancode(m_fd))
table(fish=b$fish,ritual=b$ritual2); eff(m_rit2,"x"); b$x <- b$ritual

# chadwick: kites of chiefs and kings, on the mcc tree of the global trees
g <- p[!is.na(p$lat),]
gu <- match_tips(unique(g[,c("glottocode","lineage")]),gx)
set.seed(1); gt <- maxCladeCred(lapply(gpt,graft,u=gu)); gt$edge.length[gt$edge.length<=0] <- 1e-6
g <- cbind(g,xyz(g$lat,g$lon)); g$id <- g$glottocode
h <- g[!is.na(g$jh),]; h$x <- as.vector(scale(h$jh))
H <- vcv.phylo(keep.tip(gt,unique(h$id)),corr=TRUE)
h$y <- h$ritual; m_jr <- fit("m_jh_ritual",update(m_fd,newdata=h,data2=list(V=H),seed=1,iter=4000,control=list(adapt_delta=.99)),h[,c("y","x","id","c1","c2","c3")],H,stancode(m_fd))
h$y <- h$ritual2; m_jr2 <- fit("m_jh_ritual2",update(m_fd,newdata=h,data2=list(V=H),seed=1,iter=4000,control=list(adapt_delta=.99)),h[,c("y","x","id","c1","c2","c3")],H,stancode(m_fd))
h$y <- h$contest; m_jc <- fit("m_jh_contest",update(m_fd,newdata=h,data2=list(V=H),seed=1,iter=4000,control=list(adapt_delta=.99)),h[,c("y","x","id","c1","c2","c3")],H,stancode(m_fd))
ha <- h[h$adult|h$child,]; ha$y <- ha$adult
HA <- vcv.phylo(keep.tip(gt,unique(ha$id)),corr=TRUE)
m_ja <- fit("m_jh_adult",update(m_fd,newdata=ha,data2=list(V=HA),seed=1,iter=4000,control=list(adapt_delta=.99)),ha[,c("y","x","id","c1","c2","c3")],HA,stancode(m_fd))
nrow(h); nrow(ha)
rbind(ritual=eff(m_jr,"x"),ritual_strict=eff(m_jr2,"x"),contest=eff(m_jc,"x"),adult=eff(m_ja,"x"))
aggregate(cbind(ritual,ritual2,contest)~jh,h,function(z) c(sum(z),length(z))); aggregate(adult~jh,ha,function(z) c(sum(z),length(z)))
dg <- function(m) {z <- posterior::summarise_draws(as_draws_df(m),"rhat","ess_bulk","ess_tail"); np <- nuts_params(m)
  v <- unclass(as_draws_array(m)); k <- dimnames(v)[[3]]; ph <- v[,,grep("^sd_.*Intercept",k)[1]]^2; sh <- ph/(ph+v[,,grep("^sdgp",k)[1]]^2)
  data.frame(rhat=round(max(z$rhat,na.rm=TRUE),3),worst=z$variable[which.max(z$rhat)],ess_bulk=round(min(z$ess_bulk,na.rm=TRUE)),low_ess=z$variable[which.min(z$ess_bulk)],ess_tail=round(min(z$ess_tail,na.rm=TRUE)),div=sum(np$Value[np$Parameter=="divergent__"]),
    treedepth=sum(np$Value[np$Parameter=="treedepth__"]>=10),ebfmi=round(min(rstan::get_bfmi(m$fit)),2),share_rhat=round(posterior::rhat(sh),3),share_ess=round(posterior::ess_bulk(sh)),draws=ndraws(m))}
ms <- list(m_fish=m_fish,m_fd=m_fd,m_rit=m_rit,m_rit2=m_rit2,m_jr=m_jr,m_jr2=m_jr2,m_jc=m_jc,m_ja=m_ja)
(dgt <- do.call(rbind,lapply(ms,dg))); write.csv(dgt,"diagnostics.csv")
ef <- t(sapply(list(m_fd,m_rit,m_rit2,m_jr,m_jr2,m_jc,m_ja),eff,v="x")); sh <- t(c(share(m_fish),NA)); colnames(ef) <- colnames(sh) <- c("est","lo","hi","pd")
mr <- rbind(data.frame(model=c("m_fd","m_rit","m_rit2","m_jr","m_jr2","m_jc","m_ja"),term=c("fishing dependence","ritual","ritual (strict)",rep("jurisdictional hierarchy",4)),outcome=c("fishing kite","fishing kite","fishing kite","ritual","ritual (strict)","contest","adult fliers"),ef),
  data.frame(model="m_fish",term="phylogenetic share",outcome="fishing kite",sh))
mr; write.csv(mr,"model_results.csv",row.names=FALSE)
rs <- do.call(rbind,lapply(names(res),function(f) {r <- res[[f]]; z <- sm(r); q <- r[!duplicated(r$tree),]
  data.frame(run=f,node=z$node,p=z$p[,1],p10=z$p[,2],p90=z$p[,3],age=z$age[,1],age10=z$age[,2],age90=z$age[,3],w_er=round(median(q$w_er),2),q_er=signif(median(q$q_er),2),gain=signif(median(q$gain),2),loss=signif(median(q$loss),2))}))
rs <- rs[!(rs$run=="asr_oceanic" & rs$node=="malayo_polynesian"),]
write.csv(rs,"asr_results.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"sessioninfo.txt")

# figures
w <- function(x) ifelse(x< -25,x+360,x)
wm <- map_data("world",wrap=c(-25,335))
base <- function(xl,yl) ggplot()+geom_polygon(data=wm,aes(long,lat,group=group),fill="#e6e4de",colour="white",linewidth=0.15)+
  coord_quickmap(xlim=xl,ylim=yl,expand=FALSE)+theme_void(base_size=9)+
  theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=8,colour="#2b2b2a",margin=margin(l=2,r=8)),
    legend.key.width=unit(8,"pt"),legend.margin=margin(0,0,0,0),legend.box.spacing=unit(4,"pt"),plot.margin=margin(2,2,2,2))
s$use <- ifelse(s$fishing==1,"Fishing",ifelse(s$ritual|s$myth|s$festival,"Ritual, festival, seasonal or myth",ifelse(s$present,"Other uses only","Absence reported only")))
lv <- c("Fishing","Ritual, festival, seasonal or myth","Other uses only","Absence reported only")
s$use <- factor(s$use,lv,paste0(lv," (",table(factor(s$use,lv)),")"))
f4 <- base(c(-25,335),c(-56,75))+
  geom_point(data=s[!is.na(s$lat),],aes(w(lon),lat,fill=use,shape=use),colour="#3a3a38",size=1.9,stroke=0.3,alpha=0.95)+
  scale_fill_manual(values=c("#2a78d6","#eb6834","#1baf7a","white"))+scale_shape_manual(values=c(21,24,22,21))
ggsave("fig4_map.pdf",f4,width=6.3,height=2.95); ggsave("fig4_map.png",f4,width=6.3,height=2.95,dpi=400,bg="white")
nu <- s[s$fishing==1 & s$family!="aust1307" & !is.na(s$lat),]
fm <- rbind(data.frame(lat=u$lat,lon=u$lon,state=c("Absent","Conflicting reports","Present")[match(u$fish,c(0,0.5,1))]),
  data.frame(lat=nu$lat,lon=nu$lon,state=rep("Present, non-Austronesian",nrow(nu))))
lv <- c("Present","Present, non-Austronesian","Absent","Conflicting reports")
fm <- fm[w(fm$lon)>=92 & w(fm$lon)<=256 & fm$lat>=-48 & fm$lat<=25,]
fm$state <- factor(fm$state,lv,paste0(lv," (",table(factor(fm$state,lv)),")"))
rl <- data.frame(x=c(110,150,215,168),y=c(-13,17,-30,-29),l=c("ISLAND SOUTHEAST ASIA","MICRONESIA","POLYNESIA","MELANESIA"))
f5 <- base(c(92,256),c(-48,25))+geom_text(data=rl,aes(x,y,label=l),size=2.3,colour="#8a8984",fontface="italic")+
  geom_point(data=fm,aes(w(lon),lat,fill=state,shape=state),colour="#3a3a38",size=2.1,stroke=0.3)+
  scale_fill_manual(values=c("#2a78d6","#2a78d6","white","#a8a7a2"))+scale_shape_manual(values=c(21,22,21,21))
ggsave("fig5_fishing.pdf",f5,width=6.3,height=3.25); ggsave("fig5_fishing.png",f5,width=6.3,height=3.25,dpi=400,bg="white")
set.seed(1); t1 <- graft(st,u); x <- cbind(1-u$fish,u$fish); rownames(x) <- u$glottocode; x <- x[t1$tip.label,]
an <- ancr(fitMk(t1,x,model="ER",pi="fitzjohn"))
rt <- u$ritual[match(t1$tip.label,u$glottocode)]
cn <- sapply(c(Polynesian="poly1242",Micronesian="micr1243","Central Pacific"="cent2060"),function(k) getMRCA(t1,u$glottocode[grepl(k,u$lineage)]))
t2 <- force.ultrametric(t1,method="extend")
t2$tip.label <- sapply(t1$tip.label,function(z) paste(a$society[a$glottocode==z],collapse=" / "))
pie1 <- function(x0,y0,p,r) {
  th <- seq(0,2*pi,length.out=80)
  polygon(x0+r*cos(th),y0+r*sin(th),col="#ecebe7",border=NA)
  if(p>0.004) {g <- seq(pi/2,pi/2-2*pi*min(p,1),length.out=max(4,ceiling(80*p))); polygon(c(x0,x0+r*cos(g),x0),c(y0,y0+r*sin(g),y0),col="#2a78d6",border=NA)}
  polygon(x0+r*cos(th),y0+r*sin(th),border="#3a3a38",lwd=0.5)
}
h4 <- max(nodeHeights(t2)); n4 <- Ntip(t2); lo <- 0.32
pdf(NULL); plot(t2,type="fan",open.angle=8,plot=FALSE)
pp <- base::get("last_plot.phylo",envir=.PlotPhyloEnv); dev.off()
ang <- function(p,d) {th <- atan2(p$yy[d],p$xx[d]); atan2(mean(sin(th)),mean(cos(th)))}
tips <- lapply(cn,function(k) {d <- getDescendants(t2,k); d[d<=n4]})
rot <- 90-ang(pp,tips[["Central Pacific"]])*180/pi
fan4 <- function(xl,yl,draw=TRUE) {
  plot(t2,type="fan",open.angle=8,rotate.tree=rot,cex=0.56,font=1,label.offset=lo,edge.color="#7a7975",edge.width=0.9,tip.color="#2b2b2a",x.lim=xl,y.lim=yl,no.margin=TRUE)
  p <- base::get("last_plot.phylo",envir=.PlotPhyloEnv); lw <- strwidth(t2$tip.label,cex=0.56)
  th <- atan2(p$yy[1:n4],p$xx[1:n4]); re <- h4+lo+lw
  bx <- cbind(re*cos(th),re*sin(th))
  rr <- sapply(names(cn),function(k) h4+lo+max(lw[tips[[k]]])+0.05*h4)
  rr["Central Pacific"] <- max(rr["Central Pacific"],rr["Polynesian"]+0.16*h4)
  for(k in names(cn)) {
    t <- atan2(p$yy[tips[[k]]],p$xx[tips[[k]]]); m <- atan2(mean(sin(t)),mean(cos(t)))
    dt <- (t-m+pi)%%(2*pi)-pi; g <- seq(m+min(dt)-0.02,m+max(dt)+0.02,length.out=100)
    R <- rr[k]+0.035*h4; cs <- cos(m); sn <- sin(m); tw <- strwidth(k,cex=0.72); th2 <- strheight(k,cex=0.72)
    aj <- c(if(abs(cs)<0.35) 0.5 else if(cs<0) 1 else 0,if(abs(cs)<0.35) (if(sn<0) 1 else 0) else 0.5)
    bx <- rbind(bx,cbind(rr[k]*cos(g),rr[k]*sin(g)),c(R*cs-aj[1]*tw,R*sn-aj[2]*th2),c(R*cs+(1-aj[1])*tw,R*sn+(1-aj[2])*th2))
    if(draw) {lines(rr[k]*cos(g),rr[k]*sin(g),lwd=0.9,col="#2b2b2a"); text(R*cs,R*sn,k,cex=0.72,col="#2b2b2a",adj=aj)}
  }
  if(draw) {
    ring(0.07,t2,style="ring",offset=0.12,col=ifelse(rt==1,"#eb6834","transparent"))
    for(i in seq_len(nrow(an$ace))) pie1(p$xx[n4+i],p$yy[n4+i],an$ace[i,2],0.05)
    for(i in seq_len(n4)) pie1(p$xx[i],p$yy[i],x[i,2],0.07)
    legend("topright",c("Fishing kite","No fishing kite","Ritual kite"),pch=c(21,21,NA),lty=c(0,0,1),lwd=c(1,1,3),seg.len=0.8,pt.bg=c("#2a78d6","#ecebe7",NA),col=c("#3a3a38","#3a3a38","#eb6834"),pt.cex=1.3,cex=0.7,bty="n")
  }
  apply(bx,2,range)
}
W4 <- 6.5; xl <- yl <- c(-1,1)*h4*2
for(it in 1:4) {pdf(NULL,W4,W4*diff(yl)/diff(xl)); b <- fan4(xl,yl,FALSE); dev.off(); pd <- 0.03*diff(b[,1]); xl <- b[,1]+c(-pd,pd); yl <- b[,2]+c(-pd,pd)}
H4 <- W4*diff(yl)/diff(xl)
fig6 <- function() invisible(fan4(xl,yl))
pdf("fig6_tree.pdf",W4,H4); fig6(); dev.off()
png("fig6_tree.png",W4,H4,units="in",res=400); fig6(); dev.off()
ab <- c(malayo_polynesian="PMP",oceanic="POc",micronesian="PMc",central_pacific="PCP",polynesian="PPn")
st3 <- function(r) {r$fish <- pmin(pmax(r$pf,0),1); r$kite <- pmax(r$p-r$pf,0); r$none <- pmax(1-r$p,0)
  z <- aggregate(cbind(fish,kite,none)~node,r,mean); z$age <- aggregate(age~node,r,median)$age; z}
grp <- function(l) ifelse(grepl("poly1242",l),6,ifelse(grepl("cent2060",l),5,ifelse(grepl("micr1243",l),4,ifelse(grepl("nort3195|sout3173",l),3,ifelse(grepl("ocea1241",l),2,1)))))
tl <- c("Island Southeast Asia, Palau and the Marianas","New Guinea, the Bismarcks, the Solomons and Yap","Vanuatu and New Caledonia","The Carolines, Marshalls, Kiribati and Nauru","Fiji and Rotuma","Polynesia")
s3 <- with(uj,ifelse(s2&!s0&!s1,1,ifelse(s1&!s0&!s2,2,ifelse(s0&!s1&!s2,3,4))))
tc <- unclass(table(factor(grp(uj$lineage),1:6),factor(s3,1:4))); tc
ty <- 6:1; ny <- c(4.6875,3.375,3,1.5,1); ch <- list(2,3:4,NULL,5,NULL); rw <- list(1,2:3,4,5,6); G <- 4.6
cols <- c("#2a78d6","#f0b13a","#5f5e5a","#dedcd5")
pie3 <- function(x,y,p,r) {s <- diff(par("usr")[3:4])/diff(par("usr")[1:2])*par("pin")[1]/par("pin")[2]; a0 <- pi/2
  for(k in 1:3) if(p[k]>0.002) {a <- seq(a0,a0-2*pi*p[k],length.out=max(3,ceiling(200*p[k]))); polygon(c(x,x+r*cos(a)),c(y,y+r*sin(a)*s),col=cols[k],border=NA); a0 <- a0-2*pi*p[k]}
  a <- seq(0,2*pi,length.out=100); lines(x+r*cos(a),y+r*sin(a)*s,col="white")}
side <- function(z,s) {x <- if(s<0) -z$age else G+z$age; x0 <- if(s<0) 0 else G
  for(i in 1:5) {yy <- c(ty[rw[[i]]],ny[ch[[i]]]); xx <- c(rep(x0,length(rw[[i]])),x[ch[[i]]])
    segments(x[i],min(yy,ny[i]),x[i],max(yy,ny[i]),col="#8a8984",lwd=1.2,lend=1); segments(x[i],yy,xx,yy,col="#8a8984",lwd=1.2,lend=1)}
  for(i in 1:5) {pie3(x[i],ny[i],unlist(z[i,c("fish","kite","none")]),0.19); w <- 0.6*strwidth(ab[i],cex=0.55,font=3); h <- 0.7*strheight(ab[i],cex=0.55,font=3)
    rect(x[i]-w,ny[i]+0.31-h,x[i]+w,ny[i]+0.31+h,col="white",border=NA); text(x[i],ny[i]+0.31,ab[i],cex=0.55,col="#5f5e5a",font=3)}}
fig7 <- function() {
  par(mar=c(2.6,0.2,1.3,0.2),xpd=NA); plot(NA,xlim=c(-4.4,G+4.4),ylim=c(0.55,6.45),axes=FALSE,xlab="",ylab="")
  n <- rowSums(tc); u <- (G-0.7)/max(n)
  for(i in 1:6) {x0 <- G/2-n[i]*u/2; xx <- x0+c(0,cumsum(tc[i,]))*u
    segments(c(0.05,x0+n[i]*u+0.06),ty[i],c(x0-0.06,G-0.05),ty[i],col="#cfcdc6",lwd=0.8,lty=3)
    rect(head(xx,-1),ty[i]-0.075,tail(xx,-1),ty[i]+0.075,col=cols,border="white",lwd=0.6)
    text(G/2,ty[i]+0.21,paste0(tl[i]," (",n[i],")"),cex=0.56,col="#2b2b2a")}
  side(st3(rj),-1); side(st3(rgj),1)
  axis(1,at=-(4:0),labels=4:0,cex.axis=0.6,col="#8a8984",col.axis="#5f5e5a",tcl=-0.25,mgp=c(1,0.25,0),lwd=0.8)
  axis(1,at=G+0:4,labels=0:4,cex.axis=0.6,col="#8a8984",col.axis="#5f5e5a",tcl=-0.25,mgp=c(1,0.25,0),lwd=0.8)
  mtext("Thousands of years ago",1,line=1.15,at=c(-2,G+2),cex=0.58,col="#2b2b2a")
  mtext(c("Gray et al. (2009)","Bouckaert et al. (2026)"),3,line=0.2,at=c(-2,G+2),cex=0.68,font=2,col="#2b2b2a")
  legend(G/2,0.6,c("Fishing kite","Kites without the fishing kite","No kites","Uncertain (languages only)"),fill=cols,border=NA,bty="n",cex=0.56,xjust=0.5,yjust=1,x.intersp=0.5,y.intersp=0.95,text.col="#2b2b2a")
}
pdf("fig7_states.pdf",6.5,3.5); fig7(); dev.off()
png("fig7_states.png",6.5,3.5,units="in",res=400); fig7(); dev.off()

# figures 1 to 3, public domain images
rm(list=c(ls(pattern="^m_"),"ms")); invisible(gc())
options(HTTPUserAgent="kites (research figure)")
cm <- "https://upload.wikimedia.org/wikipedia/commons/thumb/"
ia <- "https://archive.org/download/"
im <- function(url,f,b=c(0,1,0,1)) {i <- image_read(get(url,f)); w <- image_info(i)$width; h <- image_info(i)$height
  image_crop(i,geometry_area(round((b[2]-b[1])*w),round((b[4]-b[3])*h),round(b[1]*w),round(b[3]*h)))}
ink <- function(i) image_level(image_convert(i,colorspace="gray"),23.5,75.3,0.77)
gr <- function(i) image_normalize(image_convert(i,colorspace="gray"))
wipe <- function(i,b) {w <- image_info(i)$width; h <- image_info(i)$height
  for(z in b) i <- image_composite(i,image_blank(round((z[2]-z[1])*w),round((z[4]-z[3])*h),"white"),offset=paste0("+",round(z[1]*w),"+",round(z[3]*h))); i}
plate <- function(pl,rows,f,W=3000,g=30,t=74) {
  invisible(gc())
  ar <- sapply(pl,function(i) image_info(i)$width/image_info(i)$height)
  hs <- sapply(rows,function(r) (W-g*(length(r)-1))/sum(ar[r]))
  pg <- image_blank(W,round(sum(hs)+length(rows)*t+g*(length(rows)-1)),"white"); y <- 0
  for(k in seq_along(rows)) {x <- 0; y <- y+t
    for(j in rows[[k]]) {w <- round(ar[j]*hs[k]); h <- round(hs[k])
      pg <- image_composite(pg,image_border(image_resize(pl[[j]],paste0(w-4,"x",h-4,"!")),"#aaaaa8","2x2"),offset=paste0("+",round(x),"+",round(y)))
      pg <- image_annotate(pg,letters[j],size=62,weight=700,font="Helvetica",location=paste0("+",round(x),"+",round(y-t)))
      x <- x+w+g}
    y <- y+hs[k]+g}
  image_write(pg,f,depth=8,density="476x476"); rm(pg,pl); invisible(gc())
}
bal <- function(n,f,b) ink(im(paste0(ia,"essaysstudiespre00ridg/page/n",n,"_w1800.jpg"),f,b))
cul <- function(n,f,b) gr(im(paste0(ia,"cu31924023272424/page/n",n,"_w2400.jpg"),f,b))
plate(list(gr(image_rotate(im(paste0(ia,"prehistoryofavia181lauf/page/n120_w2400.jpg"),"laufer8.jpg",c(.062,.856,.127,.881)),90)),
  im(paste0(cm,"d/db/%E9%81%93%E5%85%89%E7%9A%87%E5%B8%9D%E8%A1%8C%E4%B9%90%E5%9B%BE.jpg/1920px-%E9%81%93%E5%85%89%E7%9A%87%E5%B8%9D%E8%A1%8C%E4%B9%90%E5%9B%BE.jpg"),"qing.jpg",c(.3125,.5573,0,.8547)),
  im(paste0(ia,"cu31924023272424/page/n56_w2400.jpg"),"culin.jpg",c(.1,.9005,.1075,.8384)),
  im(paste0(cm,"4/43/Kite_Flying_by_Suzuki_Harunobu_%28%E9%88%B4%E6%9C%A8_%E6%98%A5%E4%BF%A1%29.jpg/1920px-Kite_Flying_by_Suzuki_Harunobu_%28%E9%88%B4%E6%9C%A8_%E6%98%A5%E4%BF%A1%29.jpg"),"harunobu.jpg"),
  cul(63,"culin16.jpg",c(.235,.775,.215,.755)),
  im(paste0(cm,"4/43/Sikh_woodcut_on_paper_depicting_kite-flying%2C_Lahore_or_Amritsar%2C_ca.1870.jpg/1920px-Sikh_woodcut_on_paper_depicting_kite-flying%2C_Lahore_or_Amritsar%2C_ca.1870.jpg"),"sikh.jpg"),
  ink(im(paste0(ia,"prehistoryofavia181lauf/page/n122_w1800.jpg"),"bate.jpg",c(.06,.97,.195,.786))),
  im(paste0(cm,"3/3d/La_cometa.jpg/1920px-La_cometa.jpg"),"goya.jpg")),list(1:4,5:8),"fig1_plate.png")
plate(list(wipe(bal(651,"balfour586.jpg",c(.1493,.9267,.179,.4678)),list(c(.077,.092,.188,.222),c(.743,.761,.36,.386))),
  bal(662,"balfour597.jpg",c(.2379,.6562,.1247,.3443)),
  bal(671,"balfour606.jpg",c(.32,.73,.225,.49)),
  bal(650,"balfour585.jpg",c(.555,.845,.165,.52)),
  bal(652,"balfour587.jpg",c(.0724,.8385,.3069,.635)),
  bal(654,"balfour589.jpg",c(.22,.72,.22,.695)),
  bal(657,"balfour592.jpg",c(.155,.44,.38,.8)),
  bal(661,"balfour596.jpg",c(.1721,.5249,.3465,.6764)),
  bal(663,"balfour598.jpg",c(.2284,.8696,.2021,.6699))),list(1:3,4:9),"fig2_plate.png")
plate(list(gr(im(paste0(ia,"b31366648/page/n500_w2400.jpg"),"hamilton.jpg",c(.0885,.7967,.1053,.3182))),
  cul(55,"culin10.jpg",c(.235,.53,.535,.755)),
  wipe(gr(image_rotate(im(paste0(ia,"jstor-658394/page/n35_w2400.jpg"),"culin1899.jpg",c(.076,.892,.021,.98)),90)),
    list(c(.269,.284,.291,.317),c(.901,.918,.489,.522),c(.244,.258,.608,.634),c(.544,.558,.653,.679),c(.427,.445,.854,.894),c(.818,.833,.855,.888))),
  cul(62,"culin15.jpg",c(.513,.847,.147,.559)),
  im("https://iiif.wellcomecollection.org/image/V0046321/full/full/0/default.jpg","krishna.jpg",c(.1,.9,.18,.97)),
  ink(im(paste0(ia,"chepfl-lipr-AXC26/page/247_w2400.jpg"),"kestler.jpg",c(.155,.452,.16,.436)))),list(1:2,3:6),"fig3_plate.png")
