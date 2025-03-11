static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* Sample in [0,10] since (1+x)^y is mainly used for x,y > 0.
     If there was a way to distinguish between x and y, we could
     sample x in [-0.1,0.1] and y in [0,100] for example. */
  return 10 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX);
}
