      void sethash_c( int ***htab, int mm, int nchain );
      int islot( int *iv, int mm );
      int inhash_c( int *iv, int ***htab, int mm, int *label, int iconf);
      int addhash_c( int *iv, int ***htab, int mm, int label, int iconf);
      void remhash_c( int *iv, int ***htab, int mm , int label, int iconf);
      int inhash( int *iv, int **htab, int mm, int *label);
      int addhash( int *iv, int **htab, int mm, int label);
      void remhash( int *iv, int **htab, int mm , int label);
