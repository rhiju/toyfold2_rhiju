Come with calculations of `C_eff`, in as many analytical and numerical ways as you can think of, for:

1. a 2-link chain in Toyfold where each link has length of 1 and standard deviation of `sigma`, which we'll vary; assume the joint in the chain is forced to turn180°.
2. same as above for 6-link chain
3. same as above for 5-link chain (which will be under strain)
4. 2-link, 5-link, and 6-link chains allowing each joint to be either straight or bent back, with equal probability (`Delta=0` in Toyfold).

As a sanity check, I was quickly able to code up  part of this in MATLAB, and I got for 1 and 2 with `sigma=1`:

```
Number of links: 2 ==> Numerical: 0.282047  Exact: 0.282095
Number of links: 6 ==> Numerical: 0.162864  Exact: 0.162868
``` 

