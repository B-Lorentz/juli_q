function ratparse(argnum::String)
    numdenum = split(argnum, '/', keepempty=false)
    num=numdenum[1]
    if(length(numdenum)==2)
        return parse(Int64,numdenum[1])//parse(Int64,numdenum[2])
    end
    parse(Int64, num)//1
end
