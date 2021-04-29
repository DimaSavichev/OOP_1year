//
// Created by Savichev Dmitrii on 08.11.2020.
//

#ifndef BIGINTEGER_BIGINTEGER_H
#define BIGINTEGER_BIGINTEGER_H
#include <vector>
#include <string>
#include <iostream>

class BigInteger {
public:
    BigInteger() = default;
    BigInteger(const int& n);
    BigInteger(const BigInteger& number_) = default;

    int& GetSign();
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
private:
    std::vector <short int> number_;
    int sign_ = 0;
    static int const base_ = 100;

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

    Rational operator -();

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

private:
    BigInteger numerator, denominator;
    const int base_ = 100;

    void updateSign();
    void normalise();
};

//===================BigInteger=======================
BigInteger::BigInteger(const int& n) {
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
};
BigInteger& BigInteger::operator += (const BigInteger& n) {
    if (sign_ == n.sign_ || n.sign_ == 0) {
        add(n);
    } else {
        subtract(n);
    }
    return *this;
};
BigInteger& BigInteger::operator -= (const BigInteger& n) {
    if (sign_ == n.sign_ || sign_ == 0) {
        subtract(-n);
    } else {
        add(n);
    }
    return *this;
};

BigInteger& BigInteger::operator /= (const BigInteger& n) {
    if (this == &n) {
        *this = 1;
        return *this;
    }
    if (sign_ == 0 || isLessAbs(n)) {
        *this = 0;
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
        short int k = 0;
        while (n.isLessOrEqualAbs(tmp)) {
            tmp -= n;
            ++k;
        }
        ans.number_.insert(ans.number_.begin(), k);
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
        return multiply(n);
    }

};
BigInteger BigInteger::operator -() const {
    BigInteger answer = *this;
    answer.sign_ = -answer.sign_;
    return answer;
};

//======Inc&Dec======
BigInteger& BigInteger::operator ++() {
    return *this += 1;
};
BigInteger BigInteger::operator ++(int) {
    BigInteger tmp;
    ++(*this);
    return tmp;
};
BigInteger& BigInteger::operator --() {
    return *this -= 1;
};
BigInteger BigInteger::operator --(int) {
    BigInteger tmp;
    --(*this);
    return tmp;
};

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
        s = "00";
        return s;
    }
    for (long long i = number_.size() - 1; i >= 0; --i) {
        if (number_[i] / 10 == 0) {
            s += '0' + std::to_string(number_[i]);
        } else {
            s += std::to_string(number_[i]);
        }

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
BigInteger gcd(const BigInteger& first, const BigInteger& second) {
    return second ? gcd(second, first % second) : first;
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
};
Rational& Rational::operator += (const Rational& n) {
    numerator *= n.denominator;
    numerator += n.numerator * denominator;
    denominator *= n.denominator;
    normalise();
    return *this;
};
Rational& Rational::operator -= (const Rational& n) {
    numerator *= n.denominator;
    numerator -= n.numerator * denominator;
    denominator *= n.denominator;
    normalise();
    return *this;
};
Rational& Rational::operator *= (const Rational& n) {
    numerator *= n.numerator;
    denominator *= n.denominator;
    normalise();
    return *this;
};
Rational& Rational::operator /= (const Rational& n) {
    numerator *= n.denominator;
    denominator *= n.numerator;
    normalise();
    return *this;
};

Rational Rational::operator -() {
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
#endif //BIGINTEGER_BIGINTEGER_H
