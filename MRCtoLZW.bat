FOR /F %%f IN ('dir /b *.mrc') DO (
	mrc2tif -s -c lzw %%f %%~nf.tif
)