#!/usr/bin/env python
__doc__="""
$ python genTestData.py <testConfigFile> <refImpl> <sourcePath> <refPath>

Author: Stefan Dietiker <stefand@ethz.ch>
"""

import sys, os

# for debugging purposes (on python 2.x, print is not a proper function)
def myprint(s):
    print(s)
#os.system = myprint

def usage(errmsg):
    print(errmsg)
    print(__doc__)
    sys.exit(1)

def compileFileName( path, prefix, i, fileExt ):
    return os.path.join(path, "%s.%d.%s" % (prefix, i, fileExt))

def main():
    try:
        f_cfg = open(sys.argv[1])
        refImpl = os.path.abspath(sys.argv[2])
        sourcePath = os.path.abspath(sys.argv[3])
        refPath = os.path.abspath(sys.argv[4])
    except (KeyError, IOError, IndexError) as ex:
        usage(ex)
    
    lc = 0
    for line in f_cfg:
        lc += 1
        line = line.rstrip()
        if not line.startswith('#') and len(line) > 1:
            if not os.path.isfile(os.path.join(refPath,line)):
                try:
                    ar = line.split('.')
                    name, fileExt = '.'.join(ar[0:-1]), ar[-1]
                    [prefix, M, Wc, Ws, We] = name.split('-')
                    ar = prefix.split('_')
                    imgName = ar[0]
                    if len(ar) > 1:
                        scaling = ar[1]
                        # check whether scaled images already exist
                        if not os.path.isfile(os.path.join(sourcePath,prefix+'.0.'+fileExt)):
                            # call convert to scale images
                            i = 0
                            origFileName = compileFileName( sourcePath, imgName, i, fileExt )
                            while os.path.isfile(origFileName):
                                destFileName = compileFileName( sourcePath, prefix, i, fileExt )
                                os.system("convert -resize %s%% %s %s" % (
                                    scaling, origFileName, destFileName))
                                i += 1
                                origFileName = compileFileName( sourcePath, imgName, i, fileExt )
                    # call expFuse.m
                    os.system("cd %s; octave expFuse.m %s %s %s %s %s %s %s" % (
                        refImpl, sourcePath, refPath, prefix, M, Wc, Ws, We));
                        
                except ValueError:
                    print("illegal filename on line %d" % (lc,))

if __name__=="__main__":
    main()
