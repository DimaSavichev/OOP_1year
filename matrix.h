//
// Created by Savichev Dmitrii on 16.12.2020.
//

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H
#include <iostream>
#include <vector>
#define base_ 1000000000
#define base_length 9
template <unsigned N>
class Finite {
public:
    Finite(int number);
    Finite<N>& operator += (const Finite<N> second);
    Finite<N>& operator -= (const Finite<N> second);
    Finite<N>& operator *= (const Finite<N> second) ;
    Finite<N>& operator /= (const Finite<N> second);
    Finite<N>& operator ++();
    operator int() const;

    template <unsigned T>
    friend bool operator == (const Finite<T> first, const Finite<T> second);
private:
    long long number_;
    Finite<N> power(Finite<N> a, int k) const;
};

class BigInteger {
public:
    BigInteger() = default;
    BigInteger(const long long& n);
    BigInteger(const BigInteger& number_) = default;

    int& GetSign();
    int SeeSign() const;
    //========Присваивание========
    BigInteger& operator = (BigInteger integer);
    BigInteger& operator += (const BigInteger& n);
    BigInteger& operator -= (const BigInteger& n);
    BigInteger& operator /= (const BigInteger& n);
    BigInteger& operator %= (const BigInteger& n);
    BigInteger& operator *= (const BigInteger& n);

    BigInteger operator -() const ;

    //======Inc&Dec======
    BigInteger& operator ++();
    BigInteger operator ++(int);
    BigInteger& operator --();
    BigInteger operator --(int);

    //======Strings======
    std::string toString() const;
    std::string toStringSpecial() const;

    explicit operator bool() const ;

    //==========Friends===========
    friend bool operator ==(const BigInteger& first, const BigInteger& second);
    friend bool operator !=(const BigInteger& first, const BigInteger& second);
    friend bool operator >=(const BigInteger& first, const BigInteger& second);
    friend bool operator <=(const BigInteger& first, const BigInteger& second);
    friend bool operator >(const BigInteger& first, const BigInteger& second);
    friend bool operator <(const BigInteger& first, const BigInteger& second);
    friend std::ostream& operator << (std::ostream& out, const BigInteger& n);
    friend std::istream& operator >> (std::istream& in, BigInteger& n);
    friend BigInteger gcd(const BigInteger& first, const BigInteger& second);
private:
    std::vector <long long> number_;
    int sign_ = 0;
    //static long long const base_ = 1e10;

    void swap(BigInteger& integer);

    //=========Working with vectors========
    void zero();
    void resize(size_t n);
    void deleteZeros();

    //========Ariphmetic=======
    void add(const BigInteger& n);
    void subtract(const BigInteger& n);
    BigInteger& multiply(const BigInteger& n);

    //========Compare Absolute========
    bool isLessAbs(const BigInteger& n) const;
    bool isLessOrEqualAbs(const BigInteger& n) const;

    //======FasterBinaryFunctions==========
    void FastDivide(long long x);
    void FastMultiply(long long x);
};

class Rational {
public:
    Rational(const BigInteger& n): numerator(n), denominator(1) {}
    Rational(const int& n = 0): numerator(n), denominator(1) {}

    //========Присваивание========
    Rational& operator = (Rational rational);
    Rational& operator += (const Rational& n);
    Rational& operator -= (const Rational& n);
    Rational& operator *= (const Rational& n);
    Rational& operator /= (const Rational& n);

    Rational operator -() const;

    //========Strings&Doubles===========
    std::string toString() const;
    std::string asDecimal(size_t precision = 0) const;
    explicit operator double() const;

    //==========Friends===========
    friend bool operator ==(const Rational& first, const Rational& second);
    friend bool operator !=(const Rational& first, const Rational& second);
    friend bool operator >=(const Rational& first, const Rational& second);
    friend bool operator <=(const Rational& first, const Rational& second);
    friend bool operator >(const Rational& first, const Rational& second);
    friend bool operator <(const Rational& first, const Rational& second);
    friend std::ostream& operator << (std::ostream& out, const Rational& n);
    friend Rational abs(const Rational& n);
private:
    BigInteger numerator, denominator;
    //const long long base_ = 1e10;

    void updateSign();
    void normalise();
};

template <unsigned M, unsigned N, typename Field=Rational>
class Matrix {
public:
    Matrix();
    template <typename T>
    Matrix(const std::vector<std::vector<T>>& v);
    Matrix<M, N, Field>& operator += (const Matrix<M, N, Field>& second);
    Matrix<M, N, Field>& operator -= (const Matrix<M, N, Field>& second);

    std::vector<Field> getRow(int i) const;
    std::vector<Field> getColumn(int j) const;

    std::vector <Field>& operator [] (int i);
    const std::vector <Field>& operator [] (int i) const;

    Matrix<M, N, Field>& operator *= (const Matrix<N, N, Field>& second);
    Matrix<N, M, Field> transposed() const;
    Field det() const;
    Field rank() const;
    Field trace() const;

    Matrix<M, N, Field> inverted() const;
    void invert();
protected:
    std::vector<std::vector<Field>> matrix_;
};

template <unsigned N,typename Field=Rational>
using SquareMatrix = Matrix<N, N, Field>;
template <unsigned N>
Finite<N>::Finite(int number) {
    number_ = number;
    number_ %= N;
    number_ += N;
    number_ %= N;
};
template <unsigned N>
Finite<N>& Finite<N>::operator += (const Finite<N> second) {
    number_ = (number_ + second.number_) % N;
    return *this;
}
template <unsigned N>
Finite<N>& Finite<N>::operator -= (const Finite<N> second) {
    number_ = ((number_ - second.number_) % N + N) % N;
    return *this;
}
template <unsigned N>
Finite<N>& Finite<N>::operator *= (const Finite<N> second) {
    number_ = (number_ * second.number_) % N;
    return *this;
}
template <unsigned N>
Finite<N>& Finite<N>::operator /= (const Finite<N> second) {
    number_ = (number_ * power(second, N - 2)) % N;
    return *this;
}
template <unsigned N>
Finite<N>& Finite<N>::operator ++() {
    return *this += 1;
}
template <unsigned N>
Finite<N>::operator int() const{
    return number_;
}
template <unsigned N>
Finite<N> operator + (const Finite<N> first, const Finite<N> second) {
    Finite tmp = first;
    tmp += second;
    return tmp;
}
template <unsigned N>
Finite<N> operator - (const Finite<N> first, const Finite<N> second) {
    Finite tmp = first;
    tmp -= second;
    return tmp;
}
template <unsigned N>
Finite<N> operator * (const Finite<N> first, const Finite<N> second) {
    Finite tmp = first;
    tmp *= second;
    return tmp;
}
template <unsigned N>
Finite<N> operator / (const Finite<N> first, const Finite<N> second) {
    Finite tmp = first;
    tmp /= second;
    return tmp;
}
template <unsigned N>
Finite<N> Finite<N>::power(Finite<N> a, int k) const {
    Finite<N> res = 1;
    while (k) {
        if (k & 1)
            res *= a;
        a *= a;
        k >>= 1;
    }
    return res;

}
template <unsigned N>
bool operator == (const Finite<N> first, const Finite<N> second) {
    return (first.number_ - second.number_) % N == 0;
}
template <unsigned N>
bool operator != (const Finite<N> first, const Finite<N> second) {
    return !(first == second);
}
//================Matrix===============
template <unsigned M, unsigned N, typename Field>
Matrix<M, N, Field>::Matrix() {
    static_assert(M == N, "Matrix is not square");
    matrix_.resize(N);
    for (unsigned int i = 0; i < N; ++i) {
        matrix_[i].resize(N, 0);
        matrix_[i][i] = 1;
    }
}

template <unsigned M, unsigned N, typename Field>
template<typename T>
Matrix<M, N, Field>::Matrix(const std::vector<std::vector<T>>& v) {
    matrix_.resize(M);
    for (unsigned int i = 0; i < v.size(); ++i) {
        for (unsigned int j = 0; j < v[i].size(); ++j) {
            matrix_[i].push_back(static_cast<Field>(v[i][j]));
        }
    }
}

template <unsigned M, unsigned N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator += (const Matrix<M, N, Field>& second) {
    for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            matrix_[i][j] += second.matrix_[i][j];
        }
    }
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator -= (const Matrix<M, N, Field>& second) {
    for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            matrix_[i][j] -= second.matrix_[i][j];
        }
    }
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const{
    std::vector<std::vector<Field>> tmp;
    tmp.resize(N);
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
            tmp[i].push_back(matrix_[j][i]);
        }
    }
    return Matrix<N, M, Field>(tmp);
}

template <unsigned M, unsigned N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getRow(int i) const {
    return matrix_[i];
}

template <unsigned M, unsigned N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getColumn(int j) const {
    std::vector<Field> tmp;
    for (unsigned int i = 0; i < M; ++i) {
        tmp.push_back(matrix_[i][j]);
    }
    return tmp;
}

template <unsigned M, unsigned N, typename Field>
Field Matrix<M, N, Field>::trace() const {
    static_assert(M == N, "Matrix is not square");
    Field tr = 0;
    for (unsigned int i = 0; i < M; ++i) {
        tr += matrix_[i][i];
    }
    return tr;
}

template <unsigned M, unsigned N, typename Field>
std::vector <Field>& Matrix<M, N, Field>::operator [] (int i) {
    return matrix_[i];
}

template <unsigned M, unsigned N, typename Field>
const std::vector <Field>& Matrix<M, N, Field>::operator [] (int i) const {
    return matrix_[i];
}

template <unsigned M, unsigned N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator *= (const Matrix<N, N, Field>& second) {
    static_assert(M == N, "Matrix is not square");
    *this = *this * second;
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Field Matrix<M, N, Field>::det() const {
    //double precision = 1e-9;
    static_assert(M == N, "Matrix is not square");
    Matrix <N, N, Field> tmp = *this;
    Field det = 1;
    for (unsigned int i = 0; i < N; ++i) {
        unsigned int k = i;
        for (unsigned int j = i + 1; j < N; ++j) {
            if (tmp.matrix_[j][i] != 0) {
                k = j;
                break;
            }
        }
        if (tmp.matrix_[k][i] == 0) {
            return 0;
        }
        if (i != k) {
            std::swap(tmp.matrix_[i], tmp.matrix_[k]);
        } else {
            det = -det;
        }
        det *= tmp.matrix_[i][i];
        for (unsigned int j = i + 1; j < N; ++j) {
            Field coefficient = tmp.matrix_[j][i] / tmp.matrix_[i][i];
            if (tmp.matrix_[j][i] != 0) {
                for (unsigned int t = i + 1; t < N; ++t) {
                    tmp.matrix_[j][t] -= coefficient * tmp.matrix_[i][t];
                }
            }
        }

    }
    return det;
}

template <unsigned M, unsigned N, typename Field>
Field Matrix<M, N, Field>::rank() const {
//        double precision = 1e-9;
    std::vector <bool> used_lines(M);
    Matrix <M, N, Field> tmp = *this;
    Field rank = N;
    for (unsigned int i = 0; i < N; ++i) {
        unsigned int j;
        for (j = 0; j < M; ++j) {
            if (!used_lines[j] && tmp.matrix_[j][i] != 0) {
                break;
            }
        }
        if (j == M) {
            rank -= 1;
        } else {
            used_lines[j] = true;
            for (unsigned int t = i + 1; t < N; ++t) {
                Field coefficient = tmp.matrix_[j][t] / tmp.matrix_[j][i];
                for (unsigned int k = 0; k < M; ++k){
                    if (!used_lines[k] && tmp.matrix_[k][i] != 0) {
                        tmp.matrix_[k][t] -= coefficient * tmp.matrix_[k][i];
                    }
                }
            }
        }
    }
    return rank;
}

template <unsigned M, unsigned N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
//        double precision = 1e-9;
    static_assert(M == N, "Matrix is not square");
    Matrix<N, N, Field> tmp = *this;
    Matrix<N, N, Field> inverted;
    for (unsigned int i = 0; i < N; ++i) {
        unsigned int k = i;
        for (unsigned int j = i + 1; j < N; ++j) {
            if (tmp.matrix_[j][i] != 0) {
                k = j;
                break;
            }
        }
        std::swap(tmp.matrix_[i], tmp.matrix_[k]);
        std::swap(inverted.matrix_[i], inverted.matrix_[k]);
        for (unsigned int j = i + 1; j < N; ++j) {
            Field coefficient = tmp.matrix_[j][i] / tmp.matrix_[i][i];
            if (tmp.matrix_[j][i] != 0) {
                for (unsigned int t = 0; t < N; ++t) {
                    tmp.matrix_[j][t] -= coefficient * tmp.matrix_[i][t];
                    inverted.matrix_[j][t] -= coefficient * inverted.matrix_[i][t];
                }
            }
        }
    }
    for (unsigned int i = M - 1;; --i) {
        Field coefficient = tmp.matrix_[i][i];
        for (unsigned int k = 0; k < M; ++k) {
            tmp[i][k] /= coefficient;
            inverted[i][k] /= coefficient;
        }
        if (i == 0) break;
        for (unsigned int j = i - 1;; --j) {
            coefficient = tmp.matrix_[j][i];
            for (unsigned int k = 0; k < M; ++k) {
                tmp[j][k] -= coefficient * tmp.matrix_[i][k];
                inverted[j][k] -= coefficient * inverted.matrix_[i][k];
            }
            if (j == 0) break;
        }
    }
    return inverted;
}
template <unsigned M, unsigned N, typename Field>
void Matrix<M, N, Field>::invert() {
    *this = inverted();
}
template <unsigned M, unsigned N, typename Field=Rational>
Matrix<M, N, Field> operator * (const Matrix<M, N, Field>& matrix, const int a) {
    Matrix<M, N, Field> tmp = matrix;
    for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            tmp[i][j] *= a;
        }
    }
    return tmp;
}
template <unsigned M, unsigned N, typename Field = Rational>
std::ostream& operator << (std::ostream& out, Matrix<M, N, Field> matrix) {
    for (unsigned int i=0; i<M; i++) {
        for (unsigned int j=0; j<N; j++)
            out << " " << matrix[i][j];
        out << '\n';
    }
    return out;
}
template <unsigned M, unsigned N, typename Field=Rational>
Matrix<M, N, Field> operator * (const int a, const Matrix<M, N, Field>& matrix) {
    return matrix * a;
}

template <unsigned M, unsigned N, typename Field = Rational>
Matrix<M, N, Field> operator + (const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    Matrix <M, N, Field> tmp = first;
    tmp += second;
    return tmp;
}

template <unsigned M, unsigned N, typename Field = Rational>
Matrix<M, N, Field> operator - (const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    Matrix <M, N, Field> tmp = first;
    tmp -= second;
    return tmp;
}
template <unsigned M, unsigned N, unsigned K, typename Field = Rational>
Matrix<M, K, Field> operator * (const Matrix<M, N, Field>& first, const Matrix<N, K, Field>& second) {
    std::vector<std::vector <Field>> tmp;
    tmp.resize(M);
    for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < K; ++j) {
            Field x = 0;
            for (unsigned int t = 0; t < N; ++t) {
                x += first[i][t] * second[t][j];
            }
            tmp[i].push_back(x);
        }
    }
    return Matrix<M, K, Field>(tmp);
}

template <unsigned M, unsigned N, typename Field=Rational>
bool operator ==(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            if (first[i][j] != second[i][j]) {
                return false;
            }
        }
    }
    return true;
}
template <unsigned M, unsigned N, typename Field=Rational>
bool operator !=(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    return !(first == second);
}
//===================BigInteger=======================
BigInteger::BigInteger(const long long& n) {
    long long a = n;
    if (n == 0) {
        sign_ = 0;
        number_.push_back(0);
    } else if (n < 0) {
        sign_ = -1;
        a = -a;
    } else if (n > 0) {
        sign_ = 1;
    }
    while (a != 0) {
        number_.push_back(a % base_);
        a /= base_;
    }
}
int& BigInteger::GetSign() {
    return sign_;
}
//========Присваивание========
BigInteger& BigInteger::operator = (BigInteger integer) {
    swap(integer);
    return *this;
}
BigInteger& BigInteger::operator += (const BigInteger& n) {
    if (sign_ == n.sign_ || n.sign_ == 0) {
        add(n);
    } else {
        subtract(n);
    }
    return *this;
}
BigInteger& BigInteger::operator -= (const BigInteger& n) {
    if (sign_ == n.sign_ || sign_ == 0) {
        subtract(-n);
    } else {
        add(n);
    }
    return *this;
}


BigInteger& BigInteger::operator %= (const BigInteger& n) {
    BigInteger tmp = *this;
    tmp /= n;
    tmp *= n;
    *this -= tmp;
    return *this;
}
BigInteger& BigInteger::operator *= (const BigInteger& n) {
    sign_ = sign_ * n.sign_;
    if (sign_ == 0) {
        *this = 0;
        return *this;
    }
    if (this == &n) {
        return multiply(BigInteger(*this));
    } else {
        if (n.number_.size() == 1) {
            FastMultiply(n.number_[0]);
            return *this;
        } else if (number_.size() == 1) {
            long long tmp = number_[0];
            number_ = n.number_;
            FastMultiply(tmp);
            return *this;
        }
        return multiply(n);
    }

}
BigInteger BigInteger::operator -() const {
    BigInteger answer = *this;
    answer.sign_ = -answer.sign_;
    return answer;
}

//======Inc&Dec======
BigInteger& BigInteger::operator ++() {
    return *this += 1;
}
BigInteger BigInteger::operator ++(int) {
    BigInteger tmp;
    ++(*this);
    return tmp;
}
BigInteger& BigInteger::operator --() {
    return *this -= 1;
}
BigInteger BigInteger::operator --(int) {
    BigInteger tmp;
    --(*this);
    return tmp;
}

std::string BigInteger::toString() const{
    std::string s;
    if (sign_ == -1) {
        s += '-';
    }
    if (sign_ == 0) {
        s = '0';
        return s;
    }
    s += std::to_string(number_.back());
    for (long long i = number_.size() - 2; i >= 0; --i) {
        if (number_[i] / 10 == 0) {
            s += '0' + std::to_string(number_[i]);
        } else {
            s += std::to_string(number_[i]);
        }

    }
    return s;
}

std::string BigInteger::toStringSpecial() const{
    std::string s;
    if (sign_ == 0) {
        for (int i = 0; i < base_length; ++i) {
            s += '0';
        }
        return s;
    }
    for (long long i = number_.size() - 1; i >= 0; --i) {
        std::string d = std::to_string(number_[i]);
        for (unsigned int j = 0; j < base_length - d.length(); ++j) {
            s += '0';
        }
        s += d;
    }
    return s;
}

BigInteger::operator bool() const {
    return sign_ != 0;
}

void BigInteger::zero() {
    for (size_t i = 0; i < number_.size(); ++i) {
        number_[i] = 0;
    }
}
void BigInteger::swap(BigInteger& integer) {
    std::swap(number_, integer.number_);
    std::swap(sign_, integer.sign_);
}
void BigInteger::resize(size_t n) {
    while (number_.size() < n) {
        number_.push_back(0);
    }
}

void BigInteger::deleteZeros() {
    while (number_.back() == 0 && number_.size() != 1) {
        number_.pop_back();
    }
    if (number_.size() == 1 && number_.back() == 0) {
        sign_ = 0;
    }
}
void BigInteger::add(const BigInteger& n) {
    resize(n.number_.size());
    for (size_t i = 0; i < number_.size(); ++i) {
        if (i < n.number_.size()) {
            number_[i] += n.number_[i];
        } else {
            break;
        }
    }
    for (size_t i = 0; i < number_.size() - 1; ++i) {
        number_[i + 1] += number_[i] / base_;
        number_[i] %= base_;
    }
    if (number_.back() >= base_) {
        number_.back() %= base_;
        number_.push_back(1);
    }
    deleteZeros();
}

void BigInteger::subtract(const BigInteger& n) {
    sign_ = n.sign_;
    if (n == *this) {
        *this = 0;
        return;
    }
    bool is_less = n.sign_ < 0 ? n > *this : n < *this;
    resize(n.number_.size());
    number_.push_back(0);

    for (size_t i = 0; i < number_.size(); ++i) {
        if (i < n.number_.size()) {
            number_[i] = (is_less ? 1 : -1) * (number_[i] - n.number_[i]);
        }
        if (number_[i] < 0) {
            number_[i] += base_;
            number_[i + 1] += is_less ? -1 : 1;
        }

    }
    deleteZeros();
    if (is_less) {
        sign_ *= -1;
    }
}

BigInteger& BigInteger::multiply(const BigInteger& n) {
    BigInteger tmp = *this;
    zero();
    size_t k = 0;
    for (size_t j = 0; j < n.number_.size(); ++j) {
        resize(k + tmp.number_.size());
        for (size_t i = 0; i < tmp.number_.size() - 1; ++i) {
            number_[i + k] += n.number_[j] * tmp.number_[i];
            number_[i + k + 1] += number_[i + k] / base_;
            number_[i + k] %= base_;
        }
        number_[k + tmp.number_.size() - 1] += n.number_[j] * tmp.number_.back();
        if (number_[k + tmp.number_.size() - 1] >= base_) {
            if (k + tmp.number_.size() < number_.size()) {
                number_[k + tmp.number_.size()] += number_[k + tmp.number_.size() - 1] / base_;
            } else {
                number_.push_back(number_[k + tmp.number_.size() - 1] / base_);
            }
            number_[k + tmp.number_.size() - 1] %= base_;
        }
        ++k;
    }
    for (size_t i = k + tmp.number_.size() - 1; i < number_.size() - 1; ++i) {
        number_[i + 1] += number_[i] / base_;
        number_[i] %= base_;
    }
    if (number_.back() >= base_) {
        number_.push_back(number_.back() / base_);
        number_[number_.size() - 2] %= base_;
    }
    return *this;
}

bool BigInteger::isLessAbs(const BigInteger& n) const{
    if (number_.size() != n.number_.size()) {
        return number_.size() < n.number_.size();
    } else {
        for (unsigned long long i = number_.size() - 1;; --i) {
            if (number_[i] > n.number_[i]) {
                return false;
            } else if (number_[i] < n.number_[i]) {
                return true;
            }

            if (i == 0) {
                break;
            }
        }
        return false;
    }
}

bool BigInteger::isLessOrEqualAbs(const BigInteger& n) const{
    if (number_.size() != n.number_.size()) {
        return number_.size() < n.number_.size();
    } else {
        for (unsigned long long i = number_.size() - 1;; --i) {
            if (number_[i] > n.number_[i]) {
                return false;
            } else if (number_[i] < n.number_[i]) {
                return true;
            }

            if (i == 0) {
                break;
            }
        }
        return true;
    }
}

//===========Arithmetic=========
BigInteger operator +(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp += second;
    return tmp;
}
BigInteger operator -(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp -= second;
    return tmp;
}
BigInteger operator *(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp *= second;
    return tmp;
}
BigInteger operator /(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp /= second;
    return tmp;
}
BigInteger operator %(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp %= second;
    return tmp;
}
BigInteger& BigInteger::operator /= (const BigInteger& n) {
    if (this == &n) {
        *this = 1;
        return *this;
    }
    if (sign_ == 0 || isLessAbs(n)) {
        *this = 0;
        return *this;
    }
    if (n.number_.size() == 1) {
        FastDivide(n.number_[0]);
        sign_ *= n.sign_;
        return *this;
    }
    BigInteger tmp = 0;
    BigInteger ans;
    ans.sign_ = sign_ * n.sign_;
    tmp.sign_ = n.sign_;
    tmp.resize(n.number_.size());
    for (size_t i = 0; i < tmp.number_.size(); ++i) {
        tmp.number_[i] = number_[number_.size() - n.number_.size() + i];
    }
    if (tmp.isLessAbs(n)) {
        tmp.number_.insert(tmp.number_.begin(), number_[number_.size() - n.number_.size() - 1]);
    }
    long long t = number_.size() - tmp.number_.size() - 1;
    size_t p = number_.size() - tmp.number_.size();
    for (size_t i = 0; i <= p; ++i) {
        long long res = 0;
        long long left = 0;
        long long  right = base_;
        while (left <= right) {
            long long middle = (left + right) / 2;
            if ((middle * n).isLessOrEqualAbs(tmp)) {
                res = middle;
                left = middle + 1;
            } else {
                right = middle - 1;
            }
        }
        tmp -= res * n;
        ans.number_.insert(ans.number_.begin(), res);
        if (t >= 0) {
            if (tmp.sign_ == 0) {
                tmp.number_.front() = number_[t];
                if (number_[t] != 0) {
                    tmp.sign_ = n.sign_;
                }
            } else {
                tmp.number_.insert(tmp.number_.begin(), number_[t]);
                tmp.sign_ = n.sign_;
            }

            --t;
        }

    }
    *this = ans;
    return *this;
}
BigInteger gcd(const BigInteger& first, const BigInteger& second) {
    if (first == 0 || first == second) {
        return second;
    }
    if (second == 0) {
        return first;
    }
    if (first == 1 || second == 1) {
        return 1;
    }

    BigInteger g = 1;
    BigInteger u = first;
    BigInteger v = second;
    if (u.sign_ == -1) {
        u.sign_ = 1;
    }
    if (v.sign_ == -1) {
        v.sign_ = 1;
    }
    while (u.number_[0] % 2 == 0 && v.number_[0] % 2 == 0) {
        u /= 2;
        v /= 2;
        g *= 2;
    }
    while (u != 0) {
        if (u.number_[0] % 2 == 0) {
            u /= 2;
        }
        if (v.number_[0] % 2 == 0) {
            v /= 2;
        }
        if (u >= v) {
            u -= v;
        } else {
            v -= u;
        }
    }
    return g * v;
}
//========Compare========
bool operator ==(const BigInteger& first, const BigInteger& second) {
    if (first.sign_ != second.sign_) {
        return false;
    }
    if (first.number_.size() != second.number_.size()) {
        return false;
    } else {
        for (size_t i = 0; i < first.number_.size(); ++i) {
            if (first.number_[i] != second.number_[i]) {
                return false;
            }
        }
        return true;
    }
}
bool operator !=(const BigInteger& first, const BigInteger& second) {
    return !(first == second);
}
bool operator <(const BigInteger& first, const BigInteger& second) {
    if (first.sign_ < second.sign_) {
        return true;
    } else if (first.sign_ > second.sign_) {
        return false;
    } else if (first.sign_ < 0) {
        return -second < -first;
    }
    if (first.number_.size() < second.number_.size()) {
        return true;
    } else if (first.number_.size() > second.number_.size()) {
        return false;
    } else {
        for (size_t i = first.number_.size() - 1;; --i) {
            if (first.number_[i] < second.number_[i]) {
                return true;
            } else if (first.number_[i] > second.number_[i]) {
                return false;
            }
            if (i == 0) {
                break;
            }
        }
        return false;
    }
}
bool operator >(const BigInteger& first, const BigInteger& second) {
    return second < first;
}
bool operator <=(const BigInteger& first, const BigInteger& second) {
    return !(first > second);
}
bool operator >=(const BigInteger& first, const BigInteger& second) {
    return !(first < second);
}
std::ostream& operator << (std::ostream& out, const BigInteger& n) {
    out << n.toString();
    return out;
}
std::istream& operator >> (std::istream& in, BigInteger& n) {
    std::string s;
    long long start = 0;
    n.number_.clear();
    in >> s;
    if (s[0] == '-') {
        n.sign_ = -1;
        start = 1;
    } else if (s[0] == '+') {
        n.sign_ = 1;
        start = 1;
    } else {
        n.sign_ = 1;
    }

    for (long long i = s.length() - 1; i >= start; i -= 2) {
        if (i - 1 >= start) {
            n.number_.push_back((s[i - 1] - '0') * 10 + s[i] - '0');
        } else {
            n.number_.push_back(s[i] - '0');
        }

    }
    if (n.number_.back() == 0) {
        n.sign_ = 0;
    }
    return in;
}

void BigInteger::FastMultiply(long long int x) {
    long long remainder = 0;
    if (x == 0) {
        *this = 0;
        return;
    }
    if (x == 1) {
        return;
    }
    for (unsigned int i = 0; i < number_.size(); ++i) {
        number_[i] *= x;
        number_[i] += remainder;
        remainder = number_[i] / base_;
        number_[i] %= base_;
    }
    if (remainder != 0) {
        number_.push_back(remainder);
    }
    if (x < 0) {
        sign_ = -sign_;
    }
}

void BigInteger::FastDivide(long long int x) {
    if (*this == 0 || x == 1) {
        return;
    }
    std::vector<long long> result;
    long long remainder = 0;
    if (number_.back() < x) {
        result.resize(number_.size() - 1, 0);
        remainder = number_.back();
    } else {
        result.resize(number_.size(), 0);
    }
    int size = result.size();
    for (unsigned int i = size - 1;; --i) {
        result[i] = remainder * base_ + number_[i];
        remainder = result[i] % x;
        result[i] /= x;
        if (i == 0) break;
    }
    number_ = result;
}

//====================Rational===================
std::string Rational::toString() const{
    if (denominator == 1) {
        return numerator.toString();
    } else {
        return numerator.toString() + '/' + denominator.toString();
    }
}
//========Присваивание========
Rational& Rational::operator = (Rational rational) {
    numerator = rational.numerator;
    denominator = rational.denominator;
    return *this;
}
Rational& Rational::operator += (const Rational& n) {
    numerator *= n.denominator;
    numerator += n.numerator * denominator;
    denominator *= n.denominator;
    normalise();
    return *this;
}
Rational& Rational::operator -= (const Rational& n) {
    numerator *= n.denominator;
    numerator -= n.numerator * denominator;
    denominator *= n.denominator;
    normalise();
    return *this;
}
Rational& Rational::operator *= (const Rational& n) {
    numerator *= n.numerator;
    denominator *= n.denominator;
    normalise();
    return *this;
}
Rational& Rational::operator /= (const Rational& n) {
    numerator *= n.denominator;
    denominator *= n.numerator;
    normalise();
    return *this;
}

Rational Rational::operator -() const {
    Rational tmp = *this;
    tmp.numerator = -tmp.numerator;
    return tmp;
}
std::string Rational::asDecimal(size_t precision) const{
    std::string s;
    if (numerator < 0) {
        s += '-';
    }
    s += (numerator / denominator).toString();
    BigInteger tmp = numerator % denominator;
    if (precision != 0) {
        s += '.';
    }
    std::string d;
    tmp.GetSign() = 1;
    while (d.length() < precision){
        tmp *= base_;
        d += (tmp / denominator).toStringSpecial();
        tmp %= denominator;
    }
    if (d.length() > precision) {
        d.pop_back();
    }
    return s + d;
}
Rational::operator double() const {
    return stod(asDecimal(50));
}

void Rational::updateSign() {
    numerator.GetSign() *= denominator.GetSign();
    denominator.GetSign() = 1;
}
void Rational::normalise() {
    BigInteger d = gcd(numerator, denominator);
    numerator /= d;
    denominator /= d;
    updateSign();
}
//========Compare==========
bool operator ==(const Rational& first, const Rational& second) {
    return first.numerator == second.numerator && first.denominator == second.denominator;
}
bool operator !=(const Rational& first, const Rational& second) {
    return !(first == second);
}
bool operator <(const Rational& first, const Rational& second) {
    return first.numerator * second.denominator < first.denominator * second.numerator;
}
bool operator >(const Rational& first, const Rational& second) {
    return second < first;
}
bool operator <=(const Rational& first, const Rational& second) {
    return !(first > second);
}
bool operator >=(const Rational& first, const Rational& second) {
    return !(first < second);
}
//===========Arithmetic=========
Rational operator +(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp += second;
    return tmp;
}
Rational operator -(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp -= second;
    return tmp;
}
Rational operator *(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp *= second;
    return tmp;
}

Rational operator /(const Rational& first, const Rational& second) {

    Rational tmp = first;
    tmp /= second;
    return tmp;
}
std::ostream& operator << (std::ostream& out, const Rational& n) {
    out << n.toString();
    return out;
}
std::istream& operator >> (std::istream& in, Rational& n) {
    int a;
    in >> a;
    n = Rational(a);
    return in;
}
#endif //MATRIX_MATRIX_H
