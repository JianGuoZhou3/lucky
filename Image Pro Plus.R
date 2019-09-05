

# PushButton 40,35,100,28,"Image Folder",.ImgFolder # 图片文件夹
# DupImgId = IpWsDuplicate()
# ret = IpAppSelectDoc(DupImgId) #打开图片？
# ret = IpIniFileStr(GETSTRING, "IOD_Img_Path", Img_Path) #输入图片
# foundFile = Dir$(IpTrim(Img_Path) + "*."+ fExtension)
# ret = IpWsLoad(IpTrim(Img_Path)+TIFfileList(i),"tif") #载入某个图片
# Dim TIFfileList$() #定义


# Private Sub findTIFfile()


# 定义类
# Dim sName As String
# Dim ImgId As Integer

file:///
# 'count the number of TIF files in current folder, the number will be used to set list1$() dimention
# fExtension = "TIF"
# foundFile = Dir$(IpTrim(Img_Path) + "*."+ fExtension)
# TIFnumfile = 0
#While foundFile <>""
#TIFnumfile = TIFnumfile + 1
#foundFile = Dir$()
#Wend
#ReDim TIFfileList$(TIFnumfile+1) 'Redim Array dimention according to the no. of images found

#' First set the surrent sample to be current sample
# TIFfileList$(0) = sName+"."+fExtension
#foundFile = Dir$(IpTrim(Img_Path) + "*." + fExtension)
#i = 1
#While foundFile <>""
#TIFfileList$(i-1) = foundFile
#i = i+1
#foundFile = Dir$()
#Wend
#End Sub


ret = IpDocGet(GETACTDOC, 0, ImgId)










