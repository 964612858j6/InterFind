#Name: InterFind
#Author: Mingyuan Luan, Rongfeng Huang
#Maintainer: Mingyuan Luan <964612858@qq.com>

if(!require("data.table")) install.packages("data.table",update = F,ask = F)
if(!require("dplyr")) install.packages("dplyr",update = F,ask = F)
if(!require("splitstackshape")) install.packages("splitstackshape",update = F,ask = F)
if(!require("PMCMRplus")) install.packages("PMCMRplus",update = F,ask = F)
if(!require("parallel")) install.packages("parallel",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
if(!require("ggsci")) install.packages("ggsci",update = F,ask = F)
if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)

###Transform the functional gene alterations to pathway alterations
alt_trans_path=function(pathway_list, data, ncutoff=0){
	glue=function(x, sep=NULL){
		if (is.null(sep)){
			for (i in 1:length(x)){
				if (i==1){
					y=x[1]
				}
				else{
					y=paste0(y, x[i])
				}
			}
		}
		
		if (!is.null(sep)){
			for (i in 1:length(x)){
				if (i==1){
					y=x[1]
				}
				else{
					y=paste(y, x[i], sep=sep)
				}
			}
		}
		
		return(y)
	}

	path_name=names(pathway_list)
	fun=function(i){
		ord=as.numeric(as.matrix(na.omit(match(pathway_list[[i]], rownames(data)))))
		if (length(ord)>1){
			mat=data[ord,]
			
			for (j in 1:(dim(mat)[2])){
				if (j==1){
					num=length(which(mat[, j]==1))
				}
				else {
					num=c(num, length(which(mat[, j]==1)))
				}
			}
			num=ifelse(num>ncutoff, 1, 0)
			num=glue(num, sep="!!!")
			
			return(num)
		}
		else {
			return(NA)
		}
	}

	shan=sapply(1:length(pathway_list), fun)

	rowname=path_name[-which(is.na(shan))]
	shan=data.table::as.data.table(na.omit(shan))
	colnames(shan)="V1"
	shan=splitstackshape::cSplit(shan, "V1", "!!!")
	shan=as.matrix(shan)
	colnames(shan)=colnames(data)
	rownames(shan)=rowname
	
	return(t(shan))
}



###InterFind
InterFind=function(dat, category, response=NULL, returnAll=FALSE, fdr_KW=0.05, fdr_pair=0.05, cores=NULL, split_by="@@!%%!^^"){

	fdr_KW<<-fdr_KW; fdr_pair<<-fdr_pair
	split_by <<- split_by
	glue<<-function(x, sep=NULL){
		if (is.null(sep)){
			for (i in 1:length(x)){
				if (i==1){
					y=x[1]
				}
				else{
					y=paste0(y, x[i])
				}
			}
		}
	
		if (!is.null(sep)){
			for (i in 1:length(x)){
				if (i==1){
					y=x[1]
				}
				else{
					y=paste(y, x[i], sep=sep)
				}
			}
		}
		
		return(y)
	}
	
	
	if (is.null(response)){
		stop("response are needed.")
	}

###Prepare (match the samples)
	###split the cnv
	if (any(category=="cnv")){
		if (all(is.null(dim(dat)), mode(dat)=="list")){
			medi=dat[[which(category=="cnv")]]
		}
		if (!is.null(dim(dat))){
			medi=dat
		}
			shan=alist()
			shan[[1]]=medi; shan[[2]]=medi
			fun=function(x){
				rown=rownames(x); x=apply(x, 2, as.numeric); rownames(x)=rown; return(x)
			}
			shan=lapply(shan, fun)
			shan[[1]][which(shan[[1]]>0)]=1; shan[[1]][which(shan[[1]]<0)]=0
			shan[[2]][which(shan[[2]]<0)]=1; shan[[2]][which(shan[[2]]>0)]=0
		
		if (all(is.null(dim(dat)), mode(dat)=="list")){
			dat=dat[-which(category=="cnv")]
			medi=alist()
			for (i in 1:length(dat)){
				medi[[i]]=dat[[i]]
			}
			medi[[i+1]]=shan[[1]]
			medi[[i+2]]=shan[[2]]
			
			dat=medi; dat<<-dat
			medi=category[-which(category=="cnv")]; medi=c(medi, "cna", "cnd"); category=medi; category<<-category
		}
		
		save(dat, category, file="1.Rdata")
		
		if (!is.null(dim(dat))){
			dat=shan; dat<<-dat; medi=c("cna", "cnd"); category=medi; category<<-category
		}
	}
	if (all(is.null(dim(dat)), mode(dat)=="list")){
		for (i in 1:length(dat)){
			tag=data.table::as.data.table(cbind(rownames(dat[[i]]), 1:nrow(dat[[i]])))
			colnames(tag)=c("V1", category[i])
			if (i==1){
				tag_a=tag
			}
			else {
				tag_a=dplyr::left_join(tag_a, tag, by="V1")
			}
		}
		##add the response
		tag=data.table::as.data.table(cbind(rownames(response), 1:nrow(response)))
		colnames(tag)=c("V1", "response")
		tag_a=dplyr::left_join(tag_a, tag, by="V1")
		
		tag_a=na.omit(tag_a)
		
		tag_a=as.data.frame(tag_a)
		
		save(dat, category, file="2.Rdata")
		for (i in 1:length(dat)){
			nn=i+1
			dat[[i]]=dat[[i]][as.numeric(as.matrix(tag_a[,nn])), ]
		}
		response=response[as.numeric(as.matrix(tag_a[, ncol(tag_a)])), ]
		
		for (i in 1:length(dat)){
			colnames(dat[[i]])=paste(category[i], colnames(dat[[i]]), sep=split_by)
		}
		save(dat, category, file="3.Rdata")
		
		###inter
		if (length(dat)>1){
			compar=function(x){
				label=matrix(1:length(x)^2, ncol=length(x))
				colnames(label)=rownames(label)=x
				tag=label[upper.tri(label)]
				y=floor(tag/dim(label)[1])+1
				x=tag-dim(label)[1]*floor(tag/dim(label)[1])
				lis=alist()
				for (i in 1:length(x)){
				shan=label[x[i], y[i], drop=FALSE]
				lis[[i]]=c(rownames(shan), colnames(shan))
				}
				return(lis)
			}

			tag_inter=compar(1:length(dat))
		}
		
		###intra
		tag_intra=alist()
		for (i in 1:length(dat)){
			tag_intra[[i]]=c(i, i)
		}
		
		if (length(dat)>1){
			tag=c(tag_intra,tag_inter)
		}
		else{
			tag=tag_intra
		}
	}

	if (!is.null(dim(dat))){
		rowname=rownames(dat)
		dat=apply(dat, 2, as.numeric)
		rownames(dat)=rowname
		colnames(dat)=paste(category, colnames(dat), sep=split_by)
		
		sel=intersect(rownames(dat), rownames(response))
		dat=dat[match(sel, rownames(dat)), ]
		response=response[match(sel, rownames(response)), ]
		
		shan=alist()
		shan[[1]]=dat
		dat=shan; rm(shan); gc()
		
		tag=alist()
		tag[[1]]=c(1, 1)
		
		duan=TRUE
	}
	
###Run
	dat<<-dat
	response<<-response
	
	if (!is.null(cores)){
		cl.cores <<- cores
	}
	else{
		cl.cores <<- parallel::detectCores()
	}

	all=lapply(tag, function(num){
		num=as.numeric(as.matrix(num))
		x<<-dat[[num[1]]]; y<<-dat[[num[2]]]
		d1<<-dim(x)[2]
		
		cate_all<<-c("00", "11", "01", "10")
		fdr_KW<<-fdr_KW; fdr_pair<<-fdr_pair
		#x1=x; y1=y; response1=response
		#x=as.data.frame(x); y=as.data.frame(y); response=as.data.frame(response)
					
		library("parallel"); cl <<- makeCluster(cl.cores)
		clusterExport(cl, c("x", "y", "d1", "response", "cate_all", "fdr_KW", "fdr_pair", "glue", "split_by"))
		#clusterEvalQ(cl, c("data.table", "dplyr", "splitstackshape","PMCMRplus"))
		
		inf_p = parLapply(cl, 1:(dim(x)[2]*dim(y)[2]), function(w){
			if ((w/d1)%%1==0){
				dy=w/d1
				dx=d1
			}
			else {
				dy=ceiling(w/d1)
				dx=w-floor(w/d1)*d1
			}
			fun = function(dx, dy) {
				if (length(table(paste0(x[, dx], y[, dy])))==1){
					tt=paste0(x[, dx], y[, dy])
					p_a=NA
					fisher=c(NA, NA)
					p_s=NA
					me=rep(NA, 4); names(me)=c("00", "01", "10", "11")
				}
				if (length(table(paste0(x[, dx], y[, dy])))==2){
					tt=paste0(x[, dx], y[, dy])
					x1=as.numeric(as.matrix(response[which(paste0(x[, dx], y[, dy])==names(table(paste0(x[, dx], y[, dy])))[1]), 2]))
					x2=as.numeric(as.matrix(response[which(paste0(x[, dx], y[, dy])==names(table(paste0(x[, dx], y[, dy])))[2]), 2]))
					p_a=wilcox.test(x2, x1, exact=FALSE)$ p.value
					
					###fisher test
					f <- try(fisher.test(x[, dx], y[, dy]), silent = TRUE)
					if (class(f) == "try-error") {
						fisher=c(NA, NA)
					}
					else {
						fisher=c(f$p.val, f$estimate)
					}
					
					p_s=NA
					me=rep(NA, 4); names(me)=c("00", "01", "10", "11")##@@@这里要4个NA？
				}
				if (length(table(paste0(x[, dx], y[, dy])))>2){
					tt=paste0(x[, dx], y[, dy])
					p_a=kruskal.test(as.numeric(as.matrix(response[,2]))~as.factor(paste0(x[, dx], y[, dy])))$ p.value
					
					###fisher test
					f <- try(fisher.test(x[, dx], y[, dy]), silent = TRUE)
					if (class(f) == "try-error") {
						fisher=c(NA, NA)
					}
					else {
						fisher=c(f$p.val, f$estimate)
					}
					
					###Mean compare
					p_s=PMCMRplus::kwAllPairsNemenyiTest(as.numeric(as.matrix(response[,2])) ~ as.factor(paste0(x[, dx], y[, dy])))
						
					for (i in 1:length(table(p_s[[6]][,2]))){
						ww=which(as.character(as.matrix(p_s[[6]][,2]))==names(table(p_s[[6]][,2]))[i])
						ww=as.numeric(as.matrix(p_s[[6]][ww, 1]))
						if (i==1){
							me=mean(ww)
						}
						else{
							me=c(me, mean(ww))
						}
					}
					names(me)=names(table(p_s[[6]][,2]))
						
					###p_s
					p_s=p_s$ p.value
					nr=nrow(p_s)
					nam=glue(c(rownames(p_s), colnames(p_s), nr), sep=split_by)
						
					p_s=as.character(as.matrix(p_s))
						
					p_s=glue(p_s, sep=split_by)
					p_s=paste(p_s, nam, sep=split_by)
				}

				###统计分组情况
				cate=vector(length=length(cate_all))
				cate_jj=match(names(table(paste0(x[, dx], y[, dy]))), cate_all)
				for (jj in 1:length(cate_jj)){
					cate[cate_jj[jj]]=table(paste0(x[, dx], y[, dy]))[jj]
				}
				names(cate)=cate_all
					
				###合并分组与整体P值
				s=c(cate, p_a, fisher); names(s)[5:7]=c("p_KW", "p_fisher", "OR")
				
				###合并基因名称和alter种类
				nam=data.table::as.data.table(rbind(colnames(x)[dx], colnames(y)[dy]))
				colnames(nam)="V1"
				nam=as.matrix(splitstackshape::cSplit(nam, "V1", split_by))
					
				nam=c(nam[1,2], nam[2,2], nam[1,1], nam[2,1])
				
				###合并mean
				me=me[match(cate_all, names(me))]
				me=glue(me, sep=split_by)
				
				###合并所有信息，不包括分组P值
				s=c(nam, s)
				###粘贴在一起
				##定义粘贴函数
				glue=function(x, sep=NULL){
					if (is.null(sep)){
						for (i in 1:length(x)){
							if (i==1){
								y=x[1]
							}
							else{
								y=paste0(y, x[i])
							}
						}
					}
						
					if (!is.null(sep)){
						for (i in 1:length(x)){
							if (i==1){
								y=x[1]
							}
							else{
								y=paste(y, x[i], sep=sep)
							}
						}
					}
					
					return(y)
				}
				##粘贴
				s=glue(s, sep=split_by)
				
				out=alist()
				out[[1]]=s
				out[[2]]=p_s
				out[[3]]=me
				
				return(out)
			}
			return(fun(dx=dx, dy=dy))
		})
		stopCluster(cl)
		
		###fdr_adjust function
		fdr_filter=function(dat, fdr_col=NULL, criterion=NULL){
			dat=as.data.frame(dat)
			as_numeric=function(x){
				return(as.numeric(as.matrix(x)))
			}
			
			if (!is.null(fdr_col)){
				if (class(try(fdr_col%%1==0, silent = TRUE))=="try-error"){
					fdr_col=which(colnames(dat)==fdr_col)
				}
				
				p=as_numeric(dat[, fdr_col])
				p=cbind(p, 1:length(p))
				colnames(p)=c("V1", "V2")
				p=p[order(p[,1]), ]
				q=p.adjust(as_numeric(p[,1]), method="fdr")
				p[,1]=q
				p=as_numeric(p[order(as_numeric(p[, 2])), 1])
				dat[, fdr_col]=p
				if (!is.null(criterion)){
					dat=dat[which(p<criterion), ]
				}
				return(dat)
			}
			else{
				p=as_numeric(dat)
				p=cbind(p, 1:length(p))
				colnames(p)=c("V1", "V2")
				p=p[order(p[,1]), ]
				q=p.adjust(as_numeric(p[,1]), method="fdr")
				p[,1]=q
				
				if (!is.null(criterion)){
					p=p[order(as_numeric(p[, 2])), ]
					p=p[which(p[,1]<criterion), ]
					return(p)
				}
				else{
					p=as_numeric(p[order(as_numeric(p[, 2])), 1])
					return(p)
				}
				
			}
		}
		
		inf_p=unlist(inf_p)

		###split the merged outcome
		ll=1:length(inf_p)
		
		infm_m=inf_p[(0:((max(ll)/3)-1))*3+1]
		infm_ps=inf_p[(0:((max(ll)/3)-1))*3+2]
		me=inf_p[(0:((max(ll)/3)-1))*3+3]
		
		rm(inf_p); gc()
		
		##which has the sub p values
		nn=which(!is.na(infm_ps))
		
		##split infm_m
		infm_m=data.table::as.data.table(infm_m)
		colnames(infm_m)="V1"
		infm_m=splitstackshape::cSplit(infm_m, "V1", split_by)
		colnames(infm_m)=c("gene1", "gene2", "type1", "type2", "00", "11", "01", "10", "p_KW", "p_fisher", "OR")
		#fdr adjust p value of K-W test
		infm_m=fdr_filter(dat=infm_m, fdr_col="p_KW")
		#fdr adjust p value of fisher test
		infm_m=fdr_filter(dat=infm_m, fdr_col="p_fisher")
		
		colnames(infm_m)[which(colnames(infm_m)=="p_KW")]="fdr_KW"
		colnames(infm_m)[which(colnames(infm_m)=="p_fisher")]="fdr_fisher"
		Event=ifelse(infm_m$OR>1, "Co_Occurance", "Mutually_Exclusive")
		infm_m=cbind(infm_m, Event)
		
		infm_m[, which(colnames(infm_m)=="fdr_KW")]
		
		##split infm_ps
		if(length(which(is.na(infm_ps)))==length(infm_ps)){
			ps=alist()
			ps[1:length(infm_ps)]=NA
		}
		else{
			infm_ps=data.table::as.data.table(infm_ps)
			colnames(infm_ps)="V1"
			infm_ps=splitstackshape::cSplit(infm_ps, "V1", split_by)
			infm_ps=as.matrix(infm_ps)
			
			ps=alist()
			ps[1:nrow(infm_ps)]=NA
			
			for (j in 1:length(nn)){
				ll=nn[j]
				len=na.omit(infm_ps[ll,])
				len=as.numeric(as.matrix(len[length(len)]))
				part=matrix(as.numeric(as.matrix(infm_ps[ll,1:(len*len)])), nrow=len)
				
				diname=infm_ps[ll, (1:(len*2))+len^2]
				diname[which(diname==1)]="01"
				diname[which(diname=="0")]="00"
				
				rownames(part)=diname[1:len]
				colnames(part)=diname[(1:len)+len]
				
				ps[[nn[j]]]=part
			}
		}
		
		###split the me
		me=data.table::as.data.table(me); colnames(me)="V1"
		me=splitstackshape::cSplit(me, "V1", split_by)
		colnames(me)=cate_all
		
		out=alist()
		out[[1]]=infm_m
		out[[2]]=ps
		out[[3]]=me
			
		if (returnAll){
			return(out)
		}
		else{
			sel=intersect(which(infm_m$fdr_KW<fdr_KW), which(infm_m$fdr_fisher<fdr_pair))
			out[[1]]=out[[1]][sel, ]
			out[[2]]=out[[2]][sel]
			out[[3]]=out[[3]][sel,]
			return(out)
		}
	})
	
	if (all(length(all)==1, nrow(all[[1]][[1]])==0)){
		stop("No significant pair was identified. Please change the omics data or loose the cutoff.")
	}
	
	for (i in 1:length(all)){
		if (i==1){
			info=as.matrix(all[[i]][[1]])
			p_s=all[[i]][[2]]
			me=as.matrix(all[[i]][[3]])
		}
		else{
			info=rbind(info, as.matrix(all[[i]][[1]]))
			p_s=c(p_s, all[[i]][[2]])
			me=rbind(me, as.matrix(all[[i]][[3]]))
		}
	}
	
	###rm duplicated pairs
	rm_lr=function(dat, split1="%@$!", split2="?@/%", rm_duplic=TRUE, par=FALSE, cores=NULL){
		split1<<-split1
		split2<<-split2
		tag=paste(as.character(as.matrix(dat[,1])), as.character(as.matrix(dat[,2])), sep=split1)

		t1=data.table::as.data.table(names(table(tag))); colnames(t1)="V1"
		t2=splitstackshape::cSplit(t1, "V1", split1)
		t2=data.table::as.data.table(paste(as.character(as.matrix(t2[,2])), as.character(as.matrix(t2[,1])), sep=split1))

		t1<<-as.character(as.matrix(t1))
		t2<<-as.character(as.matrix(t2))

		fun<<-function(x){
			t1=as.character(as.matrix(t1))
			if (length(which(t1==as.character(as.matrix(x))))>0){
				return(paste(which(t1==as.character(as.matrix(x))),
				t1[which(t2==x)], sep=split2))
			}
			else{
				return(NA)
			}
		}
		
		if (par){
			library(parallel)
			if (!is.null(cores)){
				cl.cores <- cores
			}
			else{
				cl.cores <- detectCores()
			}
			
			cl <<- makeCluster(cl.cores)
			clusterExport(cl, c("t1", "t2", "split1", "split2", "fun"))
			
			tt=parLapply(cl, t2, fun)
			stopCluster(cl)
		}
		if (!par){
			tt=lapply(as.character(as.matrix(t2)), fun)
		}

		tt=as.character(as.matrix(tt))
		
		if (length(tt)==length(which(tt=="NA"))){
			cat("\n", "No duplicated terms", "\n", "\n")
			if (rm_duplic){
				return(dat)
			}
			else{
				return(1:nrow(dat))
			}
		}
		
		if (length(tt)!=length(which(tt=="NA"))){
			tt=data.table::as.data.table(tt); colnames(tt)="V1"
			tt=splitstackshape::cSplit(tt, "V1", split2)

			shan=1:nrow(tt)
			tt=tt[which(as.numeric(as.matrix(tt[,1]))>shan),]
			t2=t1
			t2[as.numeric(as.matrix(tt[,1]))]=as.character(as.matrix(tt[,2]))
			
			tt=data.table::as.data.table(cbind(t1, t2)); colnames(tt)=c("V1", "V2")
			tag=data.table::as.data.table(tag); colnames(tag)="V1"
			tag=dplyr::left_join(tag, tt, "V1")
			
			tag=data.table::as.data.table(as.character(as.matrix(tag[,2]))); colnames(tag)="V1"
			sel=which(!duplicated(tag))

			tag=splitstackshape::cSplit(tag, "V1", split1)
			
			if (rm_duplic){
				tag=as.data.frame(tag)
				colnames(tag)=colnames(dat)
				tag=as.matrix(tag)
				return(tag[sel, ])
			}
			else{
				return(sel)
			}
		}
		rm(split1); rm(split2); rm(t1); rm(t2); rm(fun); gc()
	}
	
	sel=cbind(
		paste(info[,1], info[,3], sep="___"),
		paste(info[,2], info[,4], sep="___")
	)

	sel=rm_lr(dat=sel, split1="%@$!", split2="?@/%", rm_duplic=FALSE)

	###outcome
	all=alist()
	all[[1]]=info[sel, ]
	all[[2]]=p_s[sel]
	all[[3]]=me[sel, ]
	
	###rm duplicated
	sel=which(paste0(as.character(as.matrix(all[[1]][,1])), as.character(as.matrix(all[[1]][,3])))!=paste0(as.character(as.matrix(all[[1]][,2])), as.character(as.matrix(all[[1]][,4]))))
	all[[1]]=all[[1]][sel, ]
	all[[2]]=all[[2]][sel]
	all[[3]]=all[[3]][sel, ]
	
	return(all)
}


analysis=function(outcome, cutoff=0.01){
	fun=function(i){
		length(which(as.numeric(as.matrix(outcome[[1]][i, 5:8]))==0))
	}

	if (nrow(outcome[[1]])>10000){
		cl <<- makeCluster(cl.cores)
		clusterExport(cl, "outcome")
		shan=as.numeric(parLapply(cl, 1:nrow(outcome[[1]]), fun))
		stopCluster(cl)
	}
	else{
		shan=sapply(1:nrow(outcome[[1]]), fun)
	}


	shan=which(shan>0)
	outcome[[1]]=outcome[[1]][-shan, ]
	outcome[[2]]=outcome[[2]][-shan]
	outcome[[3]]=outcome[[3]][-shan, ]

	###统计ps的行列数___这个不必要加进最终版本了，因为前面已经滤过了，剩下的肯定都是3了
	fun=function(i){
		return(nrow(outcome[[2]][[i]]))
	}

	if (nrow(outcome[[1]])>10000){
		cl <<- makeCluster(cl.cores)
		clusterExport(cl, "outcome")
		shan=as.numeric(parLapply(cl, 1:nrow(outcome[[1]]), fun))
		stopCluster(cl)
	}
	else{
		shan=sapply(1:nrow(outcome[[1]]), fun)
	}


	 
	classify=function(x){
		if (all(x[nrow(x), 1:ncol(x)]<cutoff)){###########1, compared with 00, 10 and 01; 2, compared with 10 and 01.
			return("yes")
		}
	}

	if (nrow(outcome[[1]])>10000){
		cl <<- makeCluster(cl.cores)
		clusterExport(cl, "cutoff")
		shan=as.character(parLapply(cl, outcome[[2]], classify))
		stopCluster(cl)
	}
	else{
		shan=as.character(lapply(outcome[[2]], classify))
	}

	table(shan)


	rspair=outcome[[1]][which(shan=="yes"), ]
	rsme=outcome[[3]][which(shan=="yes"), ]

	fun=function(i){
		return(length(which(rspair[i, 5:8]=="  0")))
	}

	fun=function(i){
		if(all(na.omit(as.numeric(as.matrix(rsme[i, c(1, 3:4)])) < as.numeric(as.matrix(rsme[i, 2]))))){###outcome[[2]] select 1:ncol(x), here we select c(1, 3:4); outcome[[2]] select 2:ncol(x), here we select c(3:4)
			return("increase")
		}
		if(all(na.omit(as.numeric(as.matrix(rsme[i, c(1, 3:4)])) > as.numeric(as.matrix(rsme[i, 2]))))){
			return("decrease")
		}
	}
	shan=sapply(1:nrow(rsme), fun)
	cat("\n", paste(length(which(shan=="increase")), "pairs can significantly increase the scores"))
	cat("\n", paste(length(which(shan=="decrease")), "pairs can significantly decrease the scores"), "\n")
	
	rspair=cbind(rspair, shan); colnames(rspair)[ncol(rspair)]="Effect"
	
	out=alist(); out[[1]]=rspair; out[[2]]=rsme
	return(out)
}



plot_analysis=function(outcome, i){
	x1=which(category==as.character(as.matrix(outcome[i, 3])))
	x2=which(category==as.character(as.matrix(outcome[i, 4])))
	col1=grep(substr(as.character(as.matrix(outcome[i, 1])), 1, 8), colnames(dat[[x1]]))
	col2=grep(substr(as.character(as.matrix(outcome[i, 2])), 1, 8), colnames(dat[[x1]]))
	
	tag=paste0(dat[[x1]][, col1], dat[[x2]][, col2])
	tag=as.data.frame(cbind(tag, response[,2]))
	tag[,1]=as.factor(tag[,1]); tag[,2]=as.numeric(as.matrix(tag[,2]))
	colnames(tag)=c("V1", "V2")
	
	source("F:/research/LAB_02/NO_00 code_depository/NO_17 compare_pair.R")
	comparsion=compar(c("00", "01", "10", "11"))
	
	p=ggplot(tag, aes(x=V1, y=V2, fill=V1))+ geom_violin(trim=TRUE, color="black") + 
		geom_boxplot(width=0.2, position=position_dodge(0.9), outlier.size = 0.5, outlier.color = "grey", color="black")+
		scale_fill_tron() +stat_compare_means(comparisons = comparsion)+ stat_compare_means(label.y = 4.5)+ guides(fill=FALSE)+
		ylab("TIDE score")
		
	p_s=PMCMRplus::kwAllPairsNemenyiTest(as.numeric(as.matrix(tag$V2)) ~ as.factor(tag$V1))
	out=alist()
	out[[1]]=p; out[[2]]=p_s
	return(out)
}

