@echo off
set PATH=%PATH%;C:\Rtools\bin;C:\Rtools\mingw_64\bin
R CMD SHLIB %*
