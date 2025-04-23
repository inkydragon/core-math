#define CORE_MATH_BIVARIATE

static inline TYPE_UNDER_TEST random_under_test_0 (void)
{
  /*
    Sample x in [-0.5,0.5] since (1+x)^y should be mainly used for x near zero.
  */
  return ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 0.5;

  /*
    Sample x in [0,10).
  */
  //  return 10 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX);
}

static inline TYPE_UNDER_TEST random_under_test_1 (void)
{
  /*
    Sample integer y in (-16,16] since (1+x)^y supposed to be used mainly for integer y.
  */
  return (rand()&0x1f) - 15.0f;

  /*
    Sample y in [0,10).
  */
  //  return 10 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX);

  /*
    Constant power y = -1.
  */
  //  return -1;
  
  /*
    Constant power y = 7/8.
  */
  //  return 7.0/8.0;
}
