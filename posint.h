#ifndef POSINT_H
#define POSINT_H

#include <iostream>
#include <vector>
#include <exception>

/* This is an exception class for the MP library. */
class MPError :public virtual std::exception {
  protected:
    const char* msg;
  public:
    MPError(const char* m = NULL) :msg(m) { }
    const char* what() const throw() 
      { return msg ? msg : "Unspecified MP error"; }
};

/* This class represents an arbitrarily large integer
 * that is at least 0. It is represented by a vector of
 * digits, starting from the least-significant digit, and
 * with each digit between 0 and B-1.
 */
class PosInt {
  private:
    // It must ALWAYS be the case that B = Bbase ^ Bpow.
    // B is really the one to be concerned about for arithmetic; 
    // Bbase just determines how the number looks for I/O operations.
    static int B;
    static int Bbase;
    static int Bpow;
   
    std::vector<int> digits;

    // Removes leading 0 digits
    void normalize();

    // Result is -1, 0, or 1 if a is <, =, or > than b,
    // up to the specified length.
    static int compareDigits (const int* a, int alen, const int* b, int blen);
    // Computes dest += x, digit-wise
    static void addArray (int* dest, const int* x, int len);
    // Computes dest -= x, digit-wise
    static void subArray (int* dest, const int* x, int len);
    // Computes dest = x * y, digit-wise
    static void mulArray 
      (int* dest, const int* x, int xlen, const int* y, int ylen);
    // Computes dest = x * y, digit-wise, using Karatsuba's method 
    // x and y must be same length
    static void fastMulArray
      (int* dest, const int* x, const int* y, int len);
    // Computes dest = dest * d, digit-wise
    static void mulDigit (int* dest, int d, int len);
    // Computes dest = dest / d, digit-wise, and returns dest % d
    static int divDigit (int* dest, int d, int len);
    // Computes division with remainder, digit-wise.
    static void divremArray 
      (int* q, int* r, const int* x, int xlen, const int* y, int ylen);

  public:
    // Computes division with remainder. After the call, we have
    // x = q*y + r, and 0 <= r < y.
    static void divrem (PosInt& q, PosInt& r, const PosInt& x, const PosInt& y);

    // Sets the base to base^pow.
    // You don't want to call this function after you've constructed
    // any PosInt objects!
    static void setBase(int base, int pow=1);

    // Default constructor. Initializes to zero
    PosInt() { }

    // Constructor from an int
    explicit PosInt (int x) { set(x); }

    // Constructor from a char array
    explicit PosInt (const char* s) { read(s); }

    // I/O routines
    void print_array(std::ostream& out) const;
    void print(std::ostream& out) const;
    void read(std::istream& in);
    void read(const char* s);

    // Sets this PosInt to the given value
    void set (int x);
    void set (const PosInt& rhs);

    // Returns this PosInt as a regular int
    int convert () const;

    // Sets this PosInt to a random number between 0 and x-1
    void rand (const PosInt& x);

    // Common comparison tests
    bool isZero() const { return digits.empty(); }
    bool isOne() const
      { return digits.size() == 1 && digits[0] == 1; }
    bool isEven() const;

    // Result is -1, 0, or 1 if this is <, =, or > than rhs.
    int compare (const PosInt& x) const;

    // this = this + x
    void add (const PosInt& x);

    // this = this - x
    void sub (const PosInt& x);

    // this = this * x
    void mul (const PosInt& x);

    // this = this * x, using Karatsuba's method
    void fastMul (const PosInt& x);

    // this = this / y
    void div (const PosInt& x)
      { PosInt temp; divrem(*this, temp, *this, x); }

    // this = this % y
    void mod (const PosInt& x)
      { PosInt temp; divrem(temp, *this, *this, x); }

    // this = this ^ x
    void pow (const PosInt& x);

    // result = a^b mod n
    void powmod (PosInt& result, const PosInt& a, const PosInt& b, const PosInt& n);

    // this = gcd(x,y)
    void gcd (const PosInt& x, const PosInt& y);

    // this = gcd(x,y) = s*x - t*y
    // NOTE THE MINUS SIGN! This is required so that both s and t are
    // always non-negative.
    void xgcd (PosInt& s, PosInt& t, const PosInt& x, const PosInt& y);

    // return true/false if this is PROBABLY prime
    bool MillerRabin () const;
};

std::ostream& operator<< (std::ostream& out, const PosInt& x);
std::istream& operator>> (std::istream& out, PosInt& x);

#endif // POSINT_H
