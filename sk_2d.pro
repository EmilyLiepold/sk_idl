
;******************************************************************************
; NAME:
;		c_sk_2D
; VERSION:
;		2.0
; PURPOSE:
;		Calculate the two-dimensional static structure factor from a sumdata file
; CATEGORY:
;		Data processing
; CALLING SEQUENCE:
;		c_sk_2D,sumdatafile,nframes=nframes
; INPUTS:
;		SUMDATAFILE:
;		a string that contents the name of the sumdatafile to be read
; OPTIONAL INPUT PARAMETERS:
;		None
; KEYWORD PARAMETERS:
;		nframes: the number of frames you want to process in the sumdata,
;					default is to process all the frames
;		diameter: particle diameter in the unit of um, default 1.58um.
;		ratio: um per pixel, default is 0.169.
;		speed: number of frames per second when the movie was grabbed. default is 30.
;		size: size of output, corresponds to resolution of the resulting plot
;		ksize, the range of k space, default is 1/ (460 pixel).
; OUTPUTS:
;		skfile that contains the calculated static structure factor.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; COMMENTS:
;		Image resolution is 1 pixel.
; MODIFICATION HISTORY:
;		Created by Bianxiao Cui on 05/16/2002.
;		Modified and commented by BC on 04/17/2004.
;		Edited by Emily Wonder on 05/16/2012 to include adjustable "size" and output options
;******************************************************************************


pro sk_2d,sumdatafile,nframes=nframes,diameter=diameter,ratio=ratio,size=size,ksize=ksize, $
skoutput=skoutput, koutput=koutput, resolution=resolution

sumdata=read_gdf(sumdatafile)
if keyword_set(nframes) then nframes=nframes else $
	nframes=max(sumdata[5,*])+1.

if keyword_set(skoutput) then skoutput=skoutput else skoutput='Sk'+strmid(sumdatafile,7)
if keyword_set(koutput) then koutput=koutput else koutput='Sk-k'

if keyword_set(size) then size=size else size=460	;make a square shaped image.

sk=fltarr(size,size)
image=intarr(size,size)
Np_total=0.
	;this part is to increase the program speed
	index=fltarr(max(sumdata[5,*])+2.)
	i=1.
	for p=0.,float(n_elements(sumdata(5,*))-1.) do begin
		if sumdata(5,p) eq i then begin
			index[i]=p
			i=i+1
		endif
	endfor
	index[i]=float(n_elements(sumdata(0,*)))
	;end of indexing sumdata file

if keyword_set(ksize) then ksize=ksize else ksize=460.
if keyword_set(resolution) then resolution=resolution else resolution=1.
nn=ksize/resolution
nn2=nn/2
T=resolution*ratio/diameter
k=(findgen(nn2)+1.)*(1./nn/T)
temp=-reverse(k)
temp=temp(1:*)
k=[0.,k,temp]*2.*3.1415926

print,nframes
;put back to 0 after this
for fr=0,nframes-1 do begin
;for some reason it's not working when fr = 1
  
	print,'frame    ',fr
	points=reform(sumdata(0:1,index[fr]:(index[fr+1]-1)))
	minx=min(points(0,*))
	miny=min(points(1,*))
	points(0,*)=points(0,*)-minx
	points(1,*)=points(1,*)-miny
	w=where(points(0,*) le size)
	points=points(*,w)
	w=where(points(1,*) le size)
	points=points(*,w)
	Np_total=Np_total+n_elements(points(0,*))
	for ptl=0,n_elements(points(0,*))-1 do begin
		x=round(points(0,*))
		y=round(points(1,*))
		image[x,y]=1.
	endfor
	ski=fft(image,1)
	sk=sk+float(ski*conj(ski))
	image=image*0
endfor
N_avg=Np_total/nframes
sk=float(sk/nframes)/N_avg

write_gdf,shift(sk,size/2,size/2),skoutput
write_gdf,k[0:(nn2-1)],koutput
end