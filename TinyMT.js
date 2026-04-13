/*========================================================================
 * TinyMT (Tiny Mersenne Twister) : A random number generator that can be
 * initialized using 3 values to create unique streams. Additionally, a
 * seed value can be used to uniformly initialize the processors.
 *
 * Usage : var tmt = new TinyMT(options) ;
 *
 * Options 
 * -------
 *      mat1  (default = 0)     : first  id
 *      mat2  (default = 0)     : second id
 *      tmat  (default = 0)     : tempering number
 *      seed  (default = 0)     : seed number for the generator
 *      linearityCheck (defaul = false) : check for linearity condition 
 *  
 *========================================================================
 */ 
class TinyMT{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CONSTRUCTOR BEGINS
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
    constructor(options={}){
        // constants .....................................................
        this._MEXP      = 127 ;
        this._SH0       = 1 ;
        this._SH1       = 10 ;
        this._SH8       = 8 ;
        this._MASK      = 0x7fffffff ;
        this._MIN_LOOP  = 8 ;
        this._PRE_LOOP  = 8 ;

        // creating the unsigned variables ...............................
        this._state = new Uint32Array(4) ;
        this._mat   = new Uint32Array(4) ;

        // reading options ...............................................
        this._mat[0] = options?.mat1 ?? 0 ;
        this._mat[1] = options?.mat2 ?? 0 ;
        this._mat[2] = options?.tmat ?? 0 ;
        this._mat[3] = options?.seed ?? 0 ;
        this._linearityCheck = options?.linearityCheck ?? false ; 

        this.spareRandom = null ;

        this.init() ;
    } // end of constructor ----------------------------------------------

/*------------------------------------------------------------------------
 * Getters and setters
 *------------------------------------------------------------------------
 */
    // getter only for constants and read-only variables ~~~~~~~~~~~~~~~~~
    get MEXP    (){ return this._MEXP       ; } 
    get SH0     (){ return this._SH0        ; }
    get SH1     (){ return this._SH1        ; }
    get SH8     (){ return this._SH8        ; }
    get MASK    (){ return this._MASK       ; }
    get MIN_LOOP(){ return this._MIN_LOOP   ; }
    get PRE_LOOP(){ return this._PRE_LOOP   ; }
    get state   (){ return this._state      ; }
    get mat     (){ return this._mat        ; }

    get mat1(){
        return this._mat[0] ;
    }
    set mat1(v){
        this.mat[0] = v ?? this.mat[0] ;
        this.init() ;
    }

    get mat2(){
        return this.mat[1] ;
    }
    set mat2(v){
        this.mat[1] = v ?? this.mat[1] ;
        this.init() ;
    }

    get tmat(){
        return this.mat[2] ;
    }
    set tmat(v){
        this.mat[2] = v ?? this.tmat ; 
        this.init() ;
    }

    get seed(){
        return this.mat[3] ;
    }
    set seed(v){
        this.mat[3] = v ?? this.seed ; 
        this.init() ;
    }


    get linearityCheck(){
        return this._linearityCheck ;
    }
    set linearityCheck(v){
        this._linearityCheck = v ;
        this.init() ;
    }

/*------------------------------------------------------------------------
 * iterate to the next state 
 *------------------------------------------------------------------------
 */
    nextState(){
        let y = this.state[3];
        let x = (this.state[0] & this.MASK)
            ^ this.state[1]
            ^ this.state[2];
        x ^= (x << this.SH0);
        y ^= (y >>> this.SH0) ^ x;
        this.state[0] = this.state[1];
        this.state[1] = this.state[2];
        this.state[2] = x ^ (y << this.SH1);
        this.state[3] = y;
        this.state[1] ^= (-(y & 1)>>>0) & this.mat1;
        this.state[2] ^= (-(y & 1)>>>0) & this.mat2;
    }
    
/*------------------------------------------------------------------------
 * initialize the generator
 *------------------------------------------------------------------------
 */
    init() {
        this.state[0] = this.seed ;
        this.state[1] = this.mat1 ;
        this.state[2] = this.mat2 ;
        this.state[3] = this.tmat ;
        for (let i = 1; i < this.MIN_LOOP; i++) {
            const  a = i & 3 ;
            const  b = (i-1) & 3 ;
            this.state[a] ^= i + Math.imul(1812433253,
                 (this.state[b]
                   ^ (this.state[b] >>> 30)));
        }

        for (let i = 0; i < this.PRE_LOOP; i++) {
            this.nextState();
        }
    }

/*------------------------------------------------------------------------
 * temper : temper the output by breaking F_2 linearity
 *------------------------------------------------------------------------
 */
    temper(){
        let t0 = new Uint32Array(1) ;
        let t1 = new Uint32Array(1) ;
        t0[0] = this.state[3];
        if (this.linearityCheck){
            t1[0] = this.state[0] ^ (this.state[2] >>> this.SH8);
        }else{
            t1[0] = this.state[0] + (this.state[2] >>> this.SH8);
        }

        t0[0] ^= t1[0] ;
        t0[0] ^= (-(t1[0] & 1)>>>0) & this.tmat;
        return t0[0] ;
    }

/*------------------------------------------------------------------------
 * randomUint32: generate an Uint32 random number
 *------------------------------------------------------------------------
 */
    randomUint32(){
        this.nextState() ;
        return this.temper() ;
    }

/*------------------------------------------------------------------------
 * randomFloat: generate a float random number between 0 and 1 
 *------------------------------------------------------------------------
 */
    randomFloat(){
        return (this.randomUint32())/4294967295. ;
    }

/*------------------------------------------------------------------------
 * random:  generate a float random number between min and max;
 *------------------------------------------------------------------------
 */
    random(min=0, max=1){
        return (this.randomFloat()*(max-min) + min) ;
    }

/*------------------------------------------------------------------------
 * randomNormal: 
 *
 * Generate random numbers that follow a Normal/Gaussian distribution
 * with a given mean and standard deviation.
 *------------------------------------------------------------------------
 */
    randomNormal(mean=0, stddev=1){
        let val, u, v, s, mul ;

        // return the spare random number if available ...................
        if (this.spareRandom !== null){
            val = this.spareRandom ;
            this.spareRandom = null ;
            return (val*stddev + mean) ;
        }

        // search for an appropriate pair ................................
        do{
            u = this.randomFloat()*2. - 1. ;
            v = this.randomFloat()*2. - 1. ;

            s = u*u + v*v ; 
        }while( s === 0 || s >= 1.) ;

        mul = Math.sqrt(-2.*Math.log(s)/s) ;

        // the pair can generate two random numbers. Save one as spare and
        // return the other one.
        this.spareRandom = v*mul ;

        return  (u*mul*stddev + mean) ;
    }

/*------------------------------------------------------------------------
 * randomNormalInRange: 
 *
 * Generate random numbers that follow a Normal distribution but 
 * are clipped to fit between min and max.
 *------------------------------------------------------------------------
 */
    randomNormalInRange(min,max,mean=0,stddev=1){
        let val ;
        do{
            val = this.randomNormal(mean,stddev) ;
        }while( val < min || val > max ) ;
        return val ;
    }

/*------------------------------------------------------------------------
 * randomLogNormal :
 *
 * Generate random numbers that follow a Log-normal distribution 
 * with a given mean and standard deviation.
 *------------------------------------------------------------------------
 */
    randomLogNormal(mean=0, stddev=1){
        return ( Math.exp( this.randomNormal(mean, stddev) ) ) ;        
    }
} /* End of TinyMT class definition ==================================== */

