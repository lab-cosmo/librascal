"Usage: vim -es < fixheader.vim ../<file_to_fix>
"Remember to set the correct header file source in the read command below
/Copyright /normal! jjd}k
read license-header.cpp.txt
wq
