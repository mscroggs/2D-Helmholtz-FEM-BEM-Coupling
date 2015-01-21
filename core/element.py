def new_linear(a,b,c):
    def func(x):
        if x<a:
            if a<b: return 0
            else: return 1
        elif x<b:
            return (x-a)/(b-a)
        elif x<c:
            return (x-c)/(b-c)
        else:
            if b<c: return 0
            else: return 1
    def derivfunc(x):
        if x<a: return 0
        elif x<b:
            return (1)/(b-a)
        elif x<c:
            return (1)/(b-c)
        else: return 0
    return func,derivfunc
