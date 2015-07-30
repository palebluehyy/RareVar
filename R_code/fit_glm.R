## Read in strand_PSEM, fit GLM model, save model
fit_poi_glm <- function(train_strand_PSEM_fnm){
	PSEM_train = read.table(strand_PSEM_fnm, as.is=T, sep='\t', header=F);
	colnames(PSEM_train) = c("chrom","pos","nalt",'depth',"transition","tripat","up_base","down_base","pre2bases", "lhlen", "lhtype",
							  "up_chdist","up_chtype","up_chlen","down_chdist","down_chtype","down_chlen","gc","hden",
							  "hrun_op","hrun_op_dist","hrun_op_len","hmer_dist","hmer_len","alt_up_down_eq");

	transit_fac_levels=c("A_A", "A_C", "A_G", "A_T", "C_A", "C_C", "C_G", "C_T", 
	"G_A", "G_C", "G_G", "G_T","T_A", "T_C", "T_G", "T_T");
	base_fac_levels=c("A","C","G","T");
	PSEM_train[,"transition"]=factor(PSEM_train[,"transition"],levels=transit_fac_levels);
	PSEM_train[,"up_base"]=factor(PSEM_train[,"up_base"],levels=base_fac_levels);
	PSEM_train[,"down_base"]=factor(PSEM_train[,"down_base"],levels=base_fac_levels);
							
	PSEM.poi.glm = glm(formula = nalt ~ transition + up_base + down_base + gc + 
    hmer_dist + hmer_len + hden + hrun_op + alt_up_down_eq + 
    offset(log(depth)), family = poisson, data = PSEM_train);
	save(PSEM.poi.glm, file=paste(train_strand_PSEM_fnm,".poi.glm.RData", sep=""));
}

predict_poi_glm <- function(glm_fnm, strand_PSEM_fnm){
	strand_PSEM = read.table(strand_PSEM_fnm, as.is=T, sep='\t', header=F);
	colnames(strand_PSEM) = c("chrom","pos","nalt",'depth',"transition","tripat","up_base","down_base","pre2bases", "lhlen", "lhtype",
							  "up_chdist","up_chtype","up_chlen","down_chdist","down_chtype","down_chlen","gc","hden",
							  "hrun_op","hrun_op_dist","hrun_op_len","hmer_dist","hmer_len","alt_up_down_eq");
	  
	transit_fac_levels=c("A_A", "A_C", "A_G", "A_T", "C_A", "C_C", "C_G", "C_T", 
	"G_A", "G_C", "G_G", "G_T","T_A", "T_C", "T_G", "T_T");
	base_fac_levels=c("A","C","G","T");
	strand_PSEM[,"transition"]=factor(strand_PSEM[,"transition"],levels=transit_fac_levels);
	strand_PSEM[,"up_base"]=factor(strand_PSEM[,"up_base"],levels=base_fac_levels);
	strand_PSEM[,"down_base"]=factor(strand_PSEM[,"down_base"],levels=base_fac_levels);
	
	load(glm_fnm);
	cat("Poisson prediction.\n");
	poi_predict=predict(PSEM.poi.glm,strand_PSEM,type="response");
	rm(PSEM.poi.glm);
	train_summary=cbind(strand_PSEM[,c("nalt","depth")],poi_predict);
	colnames(train_summary)=c("nalt","depth","poi_fit_val");
	prec_fnm = paste(strand_PSEM_fnm, ".poi.glm.prec.RData",sep="");
	result = poi_prec_loop(train_summary,0.005, 1, dim(train_summary)[1]);
	strand_PSEM = data.frame(strand_PSEM, poi_predict, result$prob_mat, result$BF_mat, stringsAsFactors=F);
	ncol = dim(strand_PSEM)[2];
	colnames(strand_PSEM)[(ncol-2) : ncol] = c("poi_prob","H1_poi_prob","poi_BF");
	save(strand_PSEM, file = prec_fnm);
}

poi_prec_loop <- function(train_summary,min_sig_af, srt_idx,end_idx){
	prob_mat=matrix(data=NA,nrow=nrow(train_summary),ncol=2);
	BF_mat=matrix(data=NA,nrow=nrow(train_summary),ncol=1);
	for(i in srt_idx:end_idx){
		mu=train_summary[i,"poi_fit_val"];
		prob_mat[i,1]=ppois(train_summary[i,"nalt"],lambda=mu);
		mu=train_summary[i,"poi_fit_val"]+train_summary[i,"depth"]*min_sig_af;
		prob_mat[i,2]=ppois(train_summary[i,"nalt"],lambda=mu);
		BF_mat[i,1]=prob_mat[i,2]/(1-prob_mat[i,1]);
	}
	colnames(prob_mat)=c("poi_prob","H1_poi_prob");
	colnames(BF_mat)=c("poi_BF");
	result = list();
	result$prob_mat = prob_mat;
	result$BF_mat = BF_mat;
	return(result);
}

FS_pvalphred_OR_totalbias<-function(vec){ ## Fisher's exact test, odds ratio, also overall strand bias: fwd/(fwd + rev)
	## format of vec: fwd_nalt, fwd_nref, rev_nalt, rev_nref
	mat=matrix(data=0,ncol=2,nrow=2);
	mat[,1]=vec[1:2];
	mat[,2]=vec[3:4];
	re=fisher.test(mat,alternative="two.sided");
	pval=log10(re$p.value)*(-10);
	OR=re$estimate;
	tb=sum(vec[1:2])/(sum(vec));
	re=c(pval,OR,tb);
	return(re);
}

filter_strand_PSEM_BF_anno <- function(fwd_strand_PSEM_fnm, rev_strand_PSEM_fnm, BF_thres, PSEM_fnm){
	## Given threshold, generate the index of loci for which BF in both strands are greater than the threshold.
	load(fwd_strand_PSEM_fnm);
	fidx = which(strand_PSEM[,'poi_BF'] >= BF_thres);
	vec_mat = cbind(strand_PSEM[,"nalt"], strand_PSEM[,"depth"] - strand_PSEM[,"nalt"]);
	rm(strand_PSEM);
	load(rev_strand_PSEM_fnm);
	ridx = which(strand_PSEM[,'poi_BF'] >= BF_thres);
	vec_mat = cbind(vec_mat, strand_PSEM[,"nalt"], strand_PSEM[,"depth"] - strand_PSEM[,"nalt"]);
	rm(strand_PSEM);	
	mat = match(fidx, ridx);
	midx = which(is.na(mat) == F);
	PSEM_idx = fidx[midx];
	
	vec_mat = vec_mat[PSEM_idx,];
	strand_re=t(apply(vec_mat,1,FS_pvalphred_OR_totalbias));
	
	PSEM = read.table(PSEM_fnm,as.is=T, sep='\t', header=F);
	colnames(PSEM) = c("chrom","pos","nalt",'depth',"AF", "transition","tripat","up_base","down_base","pre2bases", "lhlen", "lhtype",
							  "up_chdist","up_chtype","up_chlen","down_chdist","down_chtype","down_chlen","gc","hden",
							  "hrun_op","hrun_op_dist","hrun_op_len","hmer_dist","hmer_len","alt_up_down_eq");
	PSEM = PSEM[PSEM_idx,];
	
	absdiff_strandaf_percent=abs(vec_mat[,1]/rowSums(vec_mat[,c(1:2)]) -vec_mat[,3]/rowSums(vec_mat[,c(3,4)]))/PSEM[,"AF"];
	colnames(strand_re) = c("log_FS_pval","strand_OR","fwd_read_percent");
	PSEM = data.frame(PSEM, strand_re, absdiff_strandaf_percent);
	save(PSEM, file = paste(PSEM_fnm,".BF.RData", sep=""));
}


PSEM_to_vcf<-function(PSEM,qual_val,out_filenm,sample_nm){
	##########################################################################################
	## vcf input header part
	header=paste("##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	",sample_nm,sep="");
	vcf_colnames=c("chrom","pos","id","ref","alt","qual","filter","info","format","sample");
	##########################################################################################
	## vcf input body part chr_order_fnm = "/data/data/hg19_annotation/hg19.genome.ref_order.txt";
	chr_order=read.table(file=chr_order_fnm,as.is=TRUE,header=FALSE);
	chr_order=unlist(chr_order);
	chrom=PSEM[,"chrom"];
	pos=PSEM[,"pos"];
	ref=sapply(PSEM[,"transition"], function(x){strsplit(x, split="_")[[1]]});
	id=rep(".",dim(PSEM)[1]);
	alt=sapply(PSEM[,"transition"], function(x){strsplit(x, split="_")[[1]]});
	freq=PSEM[,"freq"];
	qual=rep(qual_val,dim(PSEM)[1]);
	filter=rep("PASS",dim(PSEM)[1]);
	info=paste("DP=",PSEM[,"depth"],",AF=",freq,sep="");
	format=rep("GT:DP",dim(PSEM)[1]);
	detail=paste("0/1:",PSEM[,"depth"],sep="");
	vcf=cbind(chrom,pos,id,ref,alt,qual,filter,info,format,detail);
	chrs=unique(chrom);
	mat=match(chrs,chr_order);
	mat=mat[order(mat)];
	idx=vector();
	for(i in 1:length(chrs)){
		tidx=which(vcf[,"chrom"]==chr_order[mat[i]]);
		otidx=order(as.numeric(vcf[tidx,"pos"]));
		idx=c(idx,tidx[otidx]);
	}
	new_vcf=vcf[idx,];
	out_con=file(out_filenm,"w");
	writeLines(text=header,con=out_con,sep="\n");
	close(out_con);
	write.table(new_vcf,file=out_filenm,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE);
}

extract_GATK_vcf<-function(anno_vcf_fnm,anno_vec){
	## extract annotations based on names
	## get line number to start
	cmd=paste('grep -n \"#CHROM\" ',anno_vcf_fnm,sep="");
	re=system(cmd,intern=TRUE);
	srt_line=as.numeric(strsplit(re,split=":")[[1]][1]);
	anno_vcf=read.table(anno_vcf_fnm,as.is=TRUE,sep="\t",header=FALSE,skip=srt_line);
	colnames(anno_vcf)=c("chrom","pos","ID","ref","alt","var_qual","pass","info","format","allele_depth");
	anno_mat=matrix(data=NA,ncol=length(anno_vec),nrow=dim(anno_vcf)[1]);
	for(i in 1:dim(anno_vcf)[1]){
		info=strsplit(anno_vcf[i,"info"],split=";")[[1]];
		for(a in 1:length(anno_vec)){
			idx=grep(anno_vec[a],info);
			if(length(idx)>0){
				if(anno_vec[a]=="AF"){
					anno_mat[i,a]=strsplit(info[idx],split="=")[[1]][3]
				}else{
					anno_mat[i,a]=strsplit(info[idx],split="=")[[1]][2];
				}
			}
		}
	}
	anno_mat=cbind(anno_vcf[,c("chrom","pos","ref","alt")],anno_mat);
	colnames(anno_mat)=c("chrom","pos","ref","alt",anno_vec);
	return(anno_mat);
}

rarevar_ML_features <- function(PSEM_RData_fnm, PSEM, nt){
	ref=sapply(as.character(PSEM[,"transition"]),function(x) { strsplit(x,split="_")[[1]][1] });
	alt=sapply(as.character(PSEM[,"transition"]),function(x) { strsplit(x,split="_")[[1]][2] });
	out_fnm=gsub(".RData",".inpyanno.txt",PSEM_RData_fnm);
	d=data.frame(PSEM[,c("chrom","pos")],ref,alt,stringsAsFactors=FALSE);
	colnames(d)=c("chrom","pos","ref","alt")
	write.table(d,file=out_fnm,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE);
	## add strand_wdis annotation
	strand_wdis_fnm=gsub(".inpyanno.txt",".swdis.txt",out_fnm);
	if(!file.exists(strand_wdis_fnm)){
		cmd=paste("python /rarevar/python_code/multiproc_strand_wdis_read.py -p ", nt, " ", bam_fnm," ",out_fnm," 13 S  ",strand_wdis_fnm,sep="");
		system(cmd);
	}
	swdis=read.table(strand_wdis_fnm,as.is=TRUE,sep="\t",header=TRUE);
	swdis=swdis[,c("fwd_wdis","rev_wdis")];
	## add MAPQ annotation
	MAPQ_fnm=gsub(".inpyanno.txt",".MAPQ.txt",out_fnm);
	if(!file.exists(MAPQ_fnm)){
		cmd=paste("python /rarevar/python_code/multiproc_MAPQ_by_alt_avg_diffratio.py -p ", nt, " ",bam_fnm," ",out_fnm," 13 S  ",MAPQ_fnm,sep="");
		system(cmd);
	}
	MAPQ=read.table(MAPQ_fnm,as.is=TRUE,sep="\t",header=TRUE);
	MAPQ=MAPQ[,5:dim(MAPQ)[2]]; colnames(MAPQ)=paste("MAPQ_",colnames(MAPQ),sep="");
	## add BaseQ annotation
	BaseQ_fnm=gsub(".inpyanno.txt",".BaseQ.txt",out_fnm);
	if(!file.exists(BaseQ_fnm)){
		cmd=paste("python /rarevar/python_code/multiproc_BaseQ_by_alt_avg_diffratio.py -p ", nt, " ",bam_fnm," ",out_fnm," 13 S  ",BaseQ_fnm,sep="");
		system(cmd);	
	}
	BaseQ=read.table(BaseQ_fnm,as.is=TRUE,sep="\t",header=TRUE);
	BaseQ=BaseQ[,5:dim(BaseQ)[2]]; colnames(BaseQ)=paste("BaseQ_",colnames(BaseQ),sep="");
	## add neighbour mismatches information
	nbrmm_fnm=gsub(".inpyanno.txt",".nbrmm.txt",out_fnm);
	if(!file.exists(nbrmm_fnm)){
		cmd=paste("python /rarevar/python_code/multiproc_nbr_mismatch_altpercent.py -p ", nt, " ",bam_fnm," ",out_fnm," 13 S 5 100000  ",nbrmm_fnm,sep="");
		system(cmd);
	}
	nbrmm=read.table(nbrmm_fnm,as.is=TRUE,sep="\t",header=TRUE);
	nbrmm=nbrmm[,c("mm_percent","avg_mm_perread")];
	PSEM=data.frame(PSEM,swdis,MAPQ,BaseQ,nbrmm,stringsAsFactors=FALSE);
	save(PSEM,file=PSEM_RData_fnm);

}

add_var_tag <- function(PSEM, known_var_fnm){
	## Add known variants to the PSEM. known_var_fnm is a tab-seprated file, which should contain: chrom, pos, ref, alt, af, 1st line is header.
	known_var = read.table(known_var_fnm, as.is=T, sep='\t', header=T);
	known_var_str = paste(known_var[,'chrom'], known_var[,'pos'], toupper(known_var[,'ref']), toupper(known_var[,'alt']), sep="_");
	PSEM_str = paste(PSEM[,'chrom'], PSEM[,'pos'], PSEM[,'transition'], sep="_");
	var_tag = rep("neg", dim(PSEM)[1]);
	mat = match(known_var_str, PSEM_str);
	midx = which(is.na(mat) == F);
	var_tag[mat[midx]] = "pos";
	known_var_af = rep(0, dim(PSEM)[1]);
	known_var_af[mat[midx]] = known_var[midx, "af"];
	PSEM =data.frame(PSEM, var_tag, known_var_af);
	return(PSEM);
}

hmerdistlen_anno_mat_wdis_SB_MAPQBaseQ_arff<-function(anno_mat,arff_fnm,arff_type,re_val_tag){
	cidx=which(colnames(anno_mat)=="SOR");
	if(length(cidx)>0){
		anno_mat=anno_mat[,-cidx];
	}
	idx=which(is.na(anno_mat),arr.ind=TRUE);
	if(dim(idx)[1]>0){
		cidx=unique(idx[,2]);
		for(i in 1:length(cidx)){
			ridx=which(idx[,2]==cidx[i]);
			if(is.factor(anno_mat[,cidx[i]])){
				anno_mat[,cidx[i]]=as.numeric(as.character(anno_mat[,cidx[i]]));
			}
			anno_mat[idx[ridx,1],cidx[i]]=quantile(anno_mat[,cidx[i]],probs=0.9,na.rm=TRUE);
		}
	}
	idx=vector();
	num_col=c("AF","BaseQRankSum","FS","MQ","MQRankSum","QD","ReadPosRankSum",
	"hden","hrun_op","hmer_dist","hmer_len","alt_up_down_eq","strand_OR","fwd_read_percent");	
	num_col_max=vector();
	for(c in 1:length(num_col)){
		anno_mat[,num_col[c]]=as.numeric(as.character(anno_mat[,num_col[c]]));
		tidx=which(is.infinite(anno_mat[,num_col[c]]));
		if(length(tidx)>0){
			idx=c(idx,tidx);
			num_col_max[c]=max(anno_mat[-tidx,num_col[c]]);
			anno_mat[tidx,num_col[c]]=num_col_max[c];
		}else{
			num_col_max[c]=max(anno_mat[,num_col[c]]);
		}
	}
	idx=which( (anno_mat[,"AF"]>0.95) & (rowSums(anno_mat[,c("nalt","nref")])>3000) );
	anno_mat[idx,"FS"]=0; anno_mat[idx,"strand_OR"]=1;
	
	header=paste('%
%TGen manifest file 2, pool E07set1of2
@relation \'E07set1of2\'
@attribute AlleleFreq numeric
@attribute BaseQRankSum numeric
@attribute FS numeric
@attribute MQ numeric
@attribute MQRankSum numeric
@attribute ReadPosRankSum numeric
@attribute Transition {A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_G,T_C}
@attribute up_base {"A","C","G","T","N"}
@attribute down_base {"A","C","G","T","N"}
@attribute hden numeric
@attribute hrun_op numeric
@attribute hmer_dist numeric
@attribute hmer_len numeric
@attribute alt_up_down_eq numeric
@attribute strand_OR numeric
@attribute fwd_read_percent numeric
@attribute quality_depth numeric
@attribute fwd_wdis numeric
@attribute rev_wdis numeric
@attribute absdiff_strandaf_percent numeric
@attribute MAPQ_avg_alt numeric
@attribute MAPQ_avg_ref numeric
@attribute MAPQ_alt_refdiff numeric
@attribute BaseQ_avg_alt numeric
@attribute BaseQ_avg_ref numeric
@attribute BaseQ_alt_refdiff numeric
@attribute Other_Mismatch_Percent numeric
@attribute Other_Mismatch_Perread numeric
@attribute class {yes,no}
@data',sep="");
	if(arff_type=="train"){
		idx=which(anno_mat[,"var_tag"]!="uncertain");
		anno_mat=anno_mat[idx,];
		write.table(anno_mat,file=gsub(".arff$",".noNAInf.txt",arff_fnm),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE);
		class=rep("no",dim(anno_mat)[1]);
		class[which(anno_mat[,"var_tag"]=="pos")]="yes";
		arff=cbind(anno_mat[,c("AF","BaseQRankSum","FS","MQ","MQRankSum","ReadPosRankSum","transition","up_base","down_base",
		"hden","hrun_op","hmer_dist","hmer_len","alt_up_down_eq","strand_OR","fwd_read_percent","QD","fwd_wdis","rev_wdis","absdiff_strandaf_percent",
		"MAPQ_avg_alt","MAPQ_avg_ref","MAPQ_alt_refdiff","BaseQ_avg_alt","BaseQ_avg_ref","BaseQ_alt_refdiff","mm_percent", "avg_mm_perread")],class);
	}else if(arff_type=="predict"){
		arff=cbind(anno_mat[,c("AF","BaseQRankSum","FS","MQ","MQRankSum","ReadPosRankSum","transition","up_base","down_base",
		"hden","hrun_op","hmer_dist","hmer_len","alt_up_down_eq","strand_OR","fwd_read_percent","QD","fwd_wdis","rev_wdis","absdiff_strandaf_percent",
		"MAPQ_avg_alt","MAPQ_avg_ref","MAPQ_alt_refdiff","BaseQ_avg_alt","BaseQ_avg_ref","BaseQ_alt_refdiff","mm_percent", "avg_mm_perread")],rep("?",dim(anno_mat)[1]));
	}	
	fcon=file(arff_fnm);
	writeLines(header,fcon);
	close(fcon);
	write.table(arff,file=arff_fnm,row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE,append=TRUE);
	if(re_val_tag==TRUE){
		re=list();
		re$header=header;
		re$arff=arff;
		return(re);
	}
}



training_PSEM_pipeline <- function(train_strand_PSEM_fnm, fwd_strand_PSEM_fnm, rev_strand_PSEM_fnm, BF_thres, PSEM_RData_fnm, sample_nm){
	fit_poi_glm(train_strand_PSEM_fnm); ## derive GLM model, with suffix .poi.glm.RData
	glm_fnm = paste(train_strand_PSEM_fnm,".poi.glm.RData", sep="");
	predict_poi_glm(glm_fnm, fwd_strand_PSEM_fnm); ## get GLM prediction for fwd strand, with suffix .poi.glm.prec.RData
	fwd_prec_fnm = paste(fwd_strand_PSEM_fnm, ".poi.glm.prec.RData",sep="");
	predict_poi_glm(glm_fnm, rev_strand_PSEM_fnm); ## get GLM prediction for fwd strand, with suffix .poi.glm.prec.RData
	rev_prec_fnm = paste(rev_strand_PSEM_fnm, ".poi.glm.prec.RData",sep="");
	filter_strand_PSEM_BF_anno(fwd_prec_fnm, rev_prec_fnm, BF_thres, PSEM_fnm);
	load(PSEM_RData_fnm);
	vcf_fnm = paste(gsub(".RData", "", PSEM_RData_fnm),".ML.in.vcf", sep="");
	PSEM_to_vcf(PSEM,qual_val,vcf_fnm,sample_nm);
}

run_PSEM_pipeline <- function(glm_fnm, fwd_strand_PSEM_fnm, rev_strand_PSEM_fnm, BF_thres, PSEM_RData_fnm, sample_nm){
	predict_poi_glm(glm_fnm, fwd_strand_PSEM_fnm); ## get GLM prediction for fwd strand, with suffix .poi.glm.prec.RData
	fwd_prec_fnm = paste(fwd_strand_PSEM_fnm, ".poi.glm.prec.RData",sep="");
	predict_poi_glm(glm_fnm, rev_strand_PSEM_fnm); ## get GLM prediction for fwd strand, with suffix .poi.glm.prec.RData
	rev_prec_fnm = paste(rev_strand_PSEM_fnm, ".poi.glm.prec.RData",sep="");
	filter_strand_PSEM_BF_anno(fwd_prec_fnm, rev_prec_fnm, BF_thres, PSEM_fnm);
	load(PSEM_RData_fnm);
	vcf_fnm = paste(gsub(".RData", "", PSEM_RData_fnm),".ML.in.vcf", sep="");
	PSEM_to_vcf(PSEM,qual_val,vcf_fnm,sample_nm);
}


training_ML_pipeline <- function(ref_fnm, bam_fnm, ML_vcf_fnm, bed_fnm, nt, PSEM_RData_fnm, nt){
	## generate sequencing quality annotation
	cmd=paste("/RareVar/sh/ML_GATK.sh ",ref_fnm, " ", bam_fnm," ",ML_vcf_fnm, " ", bed_fnm, " ", nt, sep="");
	system(cmd);
	ML_anno_vcf_fnm = gsub(".in.vcf", ".vcf", ML_vcf_fnm);
	## extract annotation, add to PSEM file
	anno_vec=c("AF","BaseQRankSum","Dels","FS","GC","HRun","MQ","MQRankSum","QD","ReadPosRankSum","SOR");
	anno_mat = extract_GATK_vcf(ML_anno_vcf_fnm, anno_vec);
	load(PSEM_RData_fnm);
	PSEM = data.frame(PSEM, anno_mat, stringsAsFactors=F);
	save(PSEM, file=PSEM_RData_fnm);
	rarevar_ML_features(PSEM_RData_fnm, PSEM, nt); ## add GATK and RareVar sequencing quality features
	## for training, add variant status
	load(PSEM_RData_fnm);
	PSEM = add_var_tag(PSEM, known_var_fnm);
	## generate arff file
	arff_type="train";
	re_val_tag=FALSE;
	re = hmerdistlen_anno_mat_wdis_SB_MAPQBaseQ_arff(PSEM,arff_fnm,arff_type,re_val_tag);	
	## train model with Weka random forest
}

run_ML_pipeline <- function(ref_fnm, bam_fnm, ML_vcf_fnm, bed_fnm, nt, PSEM_RData_fnm, nt){
	## generate sequencing quality annotation
	cmd=paste("/RareVar/sh/ML_GATK.sh ",ref_fnm, " ", bam_fnm," ",ML_vcf_fnm, " ", bed_fnm, " ", nt, sep="");
	system(cmd);
	ML_anno_vcf_fnm = gsub(".in.vcf", ".vcf", ML_vcf_fnm);
	## extract annotation, add to PSEM file
	anno_vec=c("BaseQRankSum","Dels","FS","GC","HRun","MQ","MQRankSum","QD","ReadPosRankSum","SOR");
	anno_mat = extract_GATK_vcf(ML_anno_vcf_fnm, anno_vec);
	load(PSEM_RData_fnm);
	PSEM = data.frame(PSEM, anno_mat, stringsAsFactors=F);
	save(PSEM, file=PSEM_RData_fnm);
	rarevar_ML_features(PSEM_RData_fnm, PSEM, nt); ## add GATK and RareVar sequencing quality features
	## generate arff file
	## run Weka random forest
}