
'My fpcov object that gives me the good dwell times should NOT be the point of maximum variance
That fpcov object, I should be able to generate by doing GRshifts plus norms
I should be able to generate the fpcov object by taking my reads, and shifting them all by (11 I think) and then resizing and getting cov.
If I was to make the fpcov object
If we do variance specific phase adjustment...
 - so the innercds adjustment actually helps slightly
 - yes AA identity does predict the A site effect
 - does taking a subset of non changing genes change my results - no.
 - is the trsums calc at least the same?
 - Why does testr have such radipaplly different numbers of okay probably multimapping
 ? What happens if I do the profile thing with only phase 0 reads
 ? Now that I have working counts by shifting I think
 -Wait fuck I think I get it. Count back from codon, 1 is 9+2
 - okay so from my codonprofiledat which works, back to the profile generation loop, then do that with overlaps
   see what reads were actually doing the overlap with, which frame of codon
 ?Ov trsums are the same?
 If yes then we check if we can get profiles that are the same
 '


#this, using the cod in subset_dwell_times, makes the correct results - high sig cor with tRNA of 29
fpcovlist<-'data/fpcovlist_oldright.rds'%>%readRDS

#but what about it's variance? Weird.
Okay, so it really doesnt identify the piont of maximum variance - it looks like the psite

goodfpcov = fpcovlist[mainsamps[1]]%>%map(~.['29'])

samp_rl_reads[[1]][['29']]%>%resize(1)%>%coverage

#these indeed resemble one another - one is missing multimappers I think
goodfpcov[[1]][[1]]%>%GRanges%>%subset(score!=0)%>%width1grs%>%resize(29)
samp_rl_reads[[1]][['29']]%>%sort

readsfromfp=goodfpcov[[1]][[1]]%>%GRanges%>%subset(score!=0)%>%width1grs%>%resize(29)

#the prestops are for sure the first nt of the prestop
cdsstarts = trspacecds%>%start%>%setNames(names(trspacecds))
cdsprestops = trspacecds%>%end%>%`-`(2+3)%>%setNames(names(trspacecds))
tr_target_pos <- cdsprestops
				# cdsprestops
GRanges(names(cdsprestops),IRanges(cdsprestops,w=6))%>%head(20)%>%exonseq[.]%>%substr(4,6)%>%{mean(is_in(.,'TGA'))>0.1}
startwinds = trspacecds%>%setstrand%>%resize(1)%>%resize(3+6)%>%resize(width(.)+6,'end')%>%
	.[!is_out_of_bounds(.)]
stopwinds = trspacecds%>%setstrand%>%resize(1,'end')%>%resize(3+6,'end')%>%resize(width(.)+6)%>%

stopwinds%>%head%>%resize(3,'center')%>%exonseq[.]

#need to run psite redo for fucntions here
readsfromfp%>%	
	trim%>%
	{.$targetpos=tr_target_pos[as.vector(seqnames(.))];.}%>%
	addphase(cdsstarts)%>%
	subset(phase==0)%>%
	{.$cor_offset=(.$targetpos+.$phase)-start(.);.}%>%
	subset(cor_offset%>%between(6,24))%>%
	.$cor_offset%>%
	table


#HOW???
scoreOv = function(gr,gr2,...){
	ov <- findOverlaps(gr, gr2, ...)
	tibble(gr=1:length(gr))%>%
		left_join(tibble(gr=ov@from,score=gr2[ov@to]$score),by='gr')%>%
		group_by(gr)%>%
		summarise(score=sum(score,na.rm=T))%>%
		.$score
}

rl=mainrls[[1]]

#YES IT"S THE SAME
ovtrsum <- lapply(mainrls,function(rl){
	scoreOv(
		innercds,
		fpcovlist[[1]][[rl]]%>%GRanges%>%subset(score!=0)%>%width1grs
	)
})%>%Reduce(`+`,.)
# quicktest(ovtrsum,trsums)

names(ovtrsum)<-names(innercds)
allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
allcodlistnz


#okay so what phase does the var offset choose? it chooses frame 1
#reduce offset to 10, increase phase to 
iphs=1
ovbasedprofscore = scoreOv(
	allcodlistnz%>%
		resize(1)%>%
		shift(45)%>%
		shift(0),
	readsfromfp%>%
		addphase(cdsstarts)%>%
		subset(phase==iphs)%>%
		resize(1)%>%
		shift(11)
)
ovbasedprofscore <- ovbasedprofscore / trsums[as.vector(seqnames(allcodlistnz))]
ovprofscore = ovbasedprofscore%>%split(cods)%>%map_dbl(mean)
ovprofscore

#YESSS okay so overlapping 
enframe(ovprofscore,'codon','ovdt')%>%left_join(
	subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-10])%>%enframe('codon','directdt')
)%>%{quicktest(.$ovdt,.$directdt)}

enframe(ovprofscore,'codon','ovdt')%>%left_join(
	codondata
)%>%{quicktest(.$ovdt,.$availability)}

subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-10])%>%enframe('codon','directdt')%>%left_join(
	codondata
)%>%{quicktest(.$directdt,.$availability)}





toptrs = seqnames(readsfromfp)%>%table%>%sort%>%tail(200)%>%names%>%intersect(names(stopwinds))

testprofreadsfromfp <- readsfromfp%>%subset(seqnames%in%toptrs)%>%addphase(cdsstarts)%>%
	subset(phase==1)%>%resize(1)%>%shift(.$offset)
testprof=testprofreadsfromfp%>%coverage%>%.[stopwinds[toptrs]]%>%as.matrix%>%colSums
testprof%>%txtplot

#HMMMMMMMMM Why is this different.
#okay let's go back to integers, no norm
seltrs = names(innercds)
rlfpcov = fpcovlist[[1]][['29']]
rlfpcov = rlfpcov[seltrs]
#rlfpcov = rlfpcov/(trsums)
# rlfpcov = rlfpcov > mean(rlfpcov)
allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
('.')
outintmats = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)
out = outintmats
for(cod in unique(cods)){
	out[[cod]] = t(out[[cod]])/trsums[as.vector(seqnames(split(allcodlistnz,cods)[[cod]]))]
}
out = out%>%map(~colMeans(t(.)))

#now why the fuck would these not match??
quicktest(subfpprofilelist[[1]][[1]][['29']][['TTT']],out[['TTT']])


################################################################################
########Okay it seems they really don't match
################################################################################
#so let's do the psite profile around codons
#
isample='E13_ribo_1'
rl='29'

bestoffsetdf<-lphaseoffsetdf%<>%filter(sample=='E13_ribo_1')%>%select(length,phase,cor_offset)
bestoffsetdf$length%<>%as.numeric

readsfromfp%<>%addphase(cdsstarts)%>%{.$length<-as.numeric(width(.));.}
readsfromfp$offset <- readsfromfp%>%mcols%>%as.data.frame%>%tibble%>%left_join(bestoffsetdf)%>%.$cor_offset

readsfromfp$seqshift<-get_probforrest_preds(
		readsfromfp%>%addphase(cdsstarts),
		seqshiftmodels[[isample]][[rl]]
	)


#
pfpcov <- readsfromfp%>%
	# shift(.$offset)%>%resize(1)%>%coverage
	shift(.$offset+.$seqshift)%>%resize(1)%>%coverage
#
pfpcov = pfpcov[seltrs]
pfpcov = pfpcov /  trsums[names(pfpcov)]
# rlfpcov = rlfpcov > mean(rlfpcov)
allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
('.')
nt=1e9
allcodlistnz2<-allcodlistnz%>%head(nt)
cods2=cods%>%head(nt)
codmat = pfpcov[allcodlistnz2]%>%split(cods2)%>%lapply(as.matrix)
out=codmat%>%map(colMeans)

# out->out_offsetonly
# out -> offsetseqshift

#now plot
plotfile<- here(paste0('plots/','psitevar_e13_1_rl29','.pdf'))
pdf(plotfile)
out%>%map_df(.id='codon',enframe,'position','count')%>%mutate(position = position-45-1)%>%
	group_by(position)%>%summarise(sd=sd(count,na.rm=T))%>%
	arrange(position)%>%
	filter(position%>%between(-18,9))%>%
	{
		qplot(data=.,x=position,y=sd)+
		theme_bw()+
		# facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('Stop/Start inferred Psite Position relative to codon')+
		# geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('blue'),linetype=2)+
		# geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of Psite occurance vs position")
	}%>%print
dev.off()
message(normalizePath(plotfile))




#
#profiles of reads if we shift by 9 - how does this looks?
testprof=readsfromfp%>%addphase(cdsstarts)%>%subset(phase==1)%>%resize(1)%>%shift(9)%>%coverage%>%.[startwinds]%>%as.matrix%>%colSums
testprof%>%txtplot

stptestprof=readsfromfp%>%addphase(cdsstarts)%>%subset(phase==1)%>%resize(1)%>%shift(9)%>%coverage%>%.[stopwinds]%>%as.matrix%>%colSums
stptestprof%>%txtplot





stopwinds%>%head%>%resize(3,'center')%>%exonseq[.]
#so I think if I take my samp reads, only the ones with phase 0, shift them by 11, I should get something
#very close to the dwell times I'm getting sucess with
psitetrsums



#same as the dwell times I have that work
samp_rl_reads[[1]][['29']]%>%subset(phase==0)


sampfpcov=fpcovlist[[1]]
rlfpcov=sampfpcov[['29']]
seltrs=names(innercds)
# trsums = sampfpcov%>%map(~.[innercds[seltrs]])%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript

rlfpcov = rlfpcov[seltrs]
rlfpcov = rlfpcov /  trsums[names(rlfpcov)]
# rlfpcov = rlfpcov > mean(rlfpcov)
allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
('.')
nt=1e9
allcodlistnz2<-allcodlistnz%>%head(nt)
cods2=cods%>%head(nt)
codmat = rlfpcov[allcodlistnz2]%>%split(cods2)%>%lapply(as.matrix)
out=codmat%>%map(colMeans)


out%>%bind_rows%>%as.matrix%>%t%>%colSds%>%txtplot((1:93)-1-45)
#okay so this gets me e.g. a 4 in the second col, now how do I get that with overlaps??
rlfpcov[allcodlistnz2%>%.[3]]
rlfpcov[allcodlistnz2%>%.[3]%>%resize(11,'start')]
readsfromfp%>%resize(1)%>%coverage(weight='score')%>%.[allcodlistnz2%>%.[3]%>%resize(11,'start')]

#
gr=allcodlistnz2%>%.[3]%>%resize(11)
gr2=readsfromfp%>%resize(1)


ovprofile <- lapply(1:93,function(i){
	scoreOv(allcodlistnz2%>%resize(1)%>%shift(i-1),readsfromfp%>%resize(1))
})



#okay so these are the values we want
subfpprofilelist[[1]][[1]][['29']]%>%map_dbl(~.[45+1-11])%>%enframe('codon','directdt')




#okay sanity at last....
(ovprofile%>%simplify2array%>%colMeans)==out[[1]]
out

foooooobar
WHY ARE THESE SO DIFFERENT

ovprofile%>%simplify2array%>%colMeans%>%`>`(0)%>%`==`(.,out[[1]]==0)
centerallcodlist <- allcodlist%>%
	resize(width(.)-3*FLANKCODS,'start')%>%
	resize(width(.)-3*FLANKCODS,'end')

#indeed taht gets us the right base pairs
rds2count <- samp_rl_reads[[1]][['29']]%>%
	# subset(phase==0)%>%
	resize(1)%>%
	subset(strand=='+')%>%
	shift(11)

trsum = innercds%>%countOverlaps(samp_rl_reads[[1]]%>%GRangesList%>%unlist)
allcoovnum <- centerallcodlist%>%countOverlaps(rds2count)
allcoovnum <- allcoovnum/trsum[as.vector(seqnames(centerallcodlist))]
allcoovnum%<>%.[is.finite(.)]
codovnum <- allcoovnum%>%split(.,names(.)%>%str_extract('\\w{3}'))%>%map_dbl(mean)

codovnum%>%enframe('codon','rdshiftov')%>%
	left_join(codondata%>%group_by(time,rep)%>%group_slice(1))%>%
	{quicktest(.$rdshiftov,.$dwell_time)}

codovnum%>%enframe('codon','rdshiftov')%>%
	left_join(codondata%>%group_by(time,rep)%>%group_slice(1))%>%
	{quicktest(.$rdshiftov,.$availability)}

#okay so that trsum is often wayyyyy too low
quicktest(log2(trsums+1),log2(trsum))
which.max((trsums+1)/(trsum+1))
testtr = 'ENSMUST00000136312'

rds2count%>%subset(seqnames%in%testtr)%>%subset()
rlfpcov[testtr]
seqlengths(rds2count)[testtr]






p1fpcov <- readsfromfp%>%	
	trim%>%
	{.$targetpos=tr_target_pos[as.vector(seqnames(.))];.}%>%
	addphase(cdsstarts)%>%
	subset(phase==0)%>%
	resize(1)
p1fpcov%<>%coverage
rlfpcov = rlfpcov > mean(rlfpcov)
allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)




