from timeit import default_timer as timer

def compare(x,y):
    if x == y:
        return True
    else:
        return False
    
time = []
i = 1
    
while i < 201:
    t0 = timer()
    x = "OIJOJFDLAKMDueoijlkajd0jewoijew0928734598u435000000000000000000000000000000000000000000000000"
    y = "OIJOJFDLAKMDueoijlkajd0jewoijew0928734598u435lkjewroijweoij98u98hjOIHJ8e9u324u92oj4o2i3jroi23"
    print(compare(x,y))
    t1 = timer()
    i += 1
    total = t1-t0
    time.append(total)

print("Avarage time: " + str(sum(time)/len(time)))


