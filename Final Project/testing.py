import os
for root, dirs, files in os.walk("..", topdown=False):
    extensions = '.csv'
    for f in files:
        ext = os.path.splitext(f)[-1]
        if ext == extensions:
            #print(os.path.join(root, f))
            print(ext)