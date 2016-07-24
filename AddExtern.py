import glob, os
os.chdir("/Users/kerblm/projects/SerpentCouplerII/src/sss")
for file in glob.glob("*.c"):
    f = open(file,'r')
    temp = f.read()
    f.close()
    if '#ifdef __cplusplus' in temp:
        print file, "already externed. "
    else:
        f = open(file, 'w')
        f.write("#ifdef __cplusplus \nextern \"C\" { \n#endif \n")
        f.write(temp)
        f.write("#ifdef __cplusplus \n} \n#endif \n")
        f.close()
        print file, "done. "
