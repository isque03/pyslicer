"""Floating-point helpers."""


class Math:
    @staticmethod
    def float_eq(a, b, eps=1e-7):
        return abs(a - b) <= eps

    # Back-compat alias used throughout legacy code paths
    floatEq = float_eq

    @staticmethod
    def remove_duplicate_intersections(mylist):
        if not mylist:
            return
        mylist.sort(key=lambda x: x.uparam)
        last = mylist[-1]
        for i in range(len(mylist) - 2, -1, -1):
            if Math.float_eq(last.uparam, mylist[i].uparam):
                del mylist[i]
            else:
                last = mylist[i]
