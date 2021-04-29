//
// Created by Savichev Dmitrii on 04.11.2020.
//
#include <iostream>
#include <cstddef>
#include <cstring>
#include <algorithm>

#ifndef STRING_STRING_H
#define STRING_STRING_H


class String {
public:
    String() = default;
    String(const char* s) {
        size_ = strlen(s);
        buffer_size_ = size_;
        str_ = new char[buffer_size_];
        memcpy(str_, s, size_);

    };
    String(const char c): String(1, c)  {};
    explicit String(size_t n, char c = '\0'): size_(n) {
        buffer_size_ = size_;
        str_ = new char[buffer_size_];
        memset(str_, c, size_);
    };
    String(const String& s): String(s.size_, '\0') {
        memcpy(str_, s.str_, size_);
    };
    ~String() {
        delete[] str_;
    }

    size_t length() const{
        return size_;
    };
    void push_back(const char& c) {
        if (size_ + 1 <= buffer_size_) {
            str_[size_] = c;
            ++size_;
        } else {
            if (buffer_size_ == 0) {
                buffer_size_ = 1;
            } else {
                buffer_size_ *= 2;
            }
            char* new_str = new char[buffer_size_];
            for (size_t i = 0; i < size_; ++i) {
                new_str[i] = str_[i];
            }
            new_str[size_] = c;
            ++size_;
            delete[] str_;
            str_ = new_str;
        }
    };
    void pop_back() {
        if (size_ - 1 >= buffer_size_ / 4) {
            --size_;
        } else {
            buffer_size_ /= 2;

            char* new_str = new char[buffer_size_];
            --size_;
            for (size_t i = 0; i < size_; ++i) {
                new_str[i] = str_[i];
            }
            delete[] str_;
            str_ = new_str;
        }
    };
    char& front() {
        return str_[0];
    };
    char& back() {
        return str_[size_ - 1];
    };
    const char& front() const{
        return str_[0];
    };
    const char& back() const {
        return str_[size_ - 1];
    };
    size_t find(const String& s) const {
        if (size_ >= s.size_) {
            for (size_t i = 0; i < size_ - s.size_ + 1; ++i) {
                bool found = true;
                for (size_t j = 0; j < s.size_; ++j) {
                    if (str_[i + j] != s.str_[j]) {
                        found = false;
                        break;
                    }
                }
                if (found) {
                    return i;
                }
            }
        }

        return size_;
    };
    size_t rfind(const String& s) const {
        for (size_t i = size_ - 1; i >= s.size_ - 1; --i) {
            bool found = true;
            for (size_t j = s.size_ - 1;; --j) {
                if (str_[i - (s.size_ - j - 1)] != s.str_[j]) {
                    found = false;
                    break;
                }
                if (j == 0) {
                    break;
                }
            }
            if (found) {
                return i - s.size_ + 1;
            }
        }
        return size_;
    };

    String substr(size_t start, size_t count) const {
        String s = String(count);
        memcpy(s.str_, str_ + start, count);
        return s;
    };
    bool empty() const {
        return size_ == 0;
    };
    void clear() {
        delete[] str_;
        size_ = 0;
        buffer_size_ = 0;
        str_ = nullptr;
    };

    char& operator [](size_t index) {
        return str_[index];
    }
    const char& operator [](size_t index) const {
        return str_[index];
    }
    String& operator =(String s) {
        swap(s);
        return *this;
    };
    String& operator +=(const String& s) {
        if (buffer_size_ >= size_ + s.size_) {
            memcpy(str_ + size_, s.str_, s.size_);
            size_ += s.size_;
        } else {
            char* new_str = new char[2 * (size_ + s.size_)];
            memcpy(new_str, str_, size_);
            memcpy(new_str + size_, s.str_, s.size_);
            delete[] str_;
            str_ = new_str;
            size_ += s.size_;
            buffer_size_ = 2 * size_;
        }
        return *this;
    };

private:
    void swap(String& s) {
        std::swap(size_, s.size_);
        std::swap(buffer_size_, s.buffer_size_);
        std::swap(str_, s.str_);
    }

    size_t size_ = 0;
    size_t buffer_size_ = 0;
    char* str_ = nullptr;

};
std::ostream& operator << (std::ostream& out, const String& s) {
    for (size_t i = 0; i < s.length(); ++i) {
        out << s[i];
    }
    return out;
}
std::istream& operator >> (std::istream& in, String& s) {
    char c;
    s.clear();
    while (in >> std::noskipws >> c && !isspace(c)) {
        s.push_back(c);
    }
    return in;
}
String operator +(const String& first, const String& second) {
    String sum = first;
    return sum += second;
}
bool operator ==(const String& first, const String& second) {
    if (first.length() != second.length()) {
        return false;
    } else {
        for (size_t i = 0; i < first.length(); ++i) {
            if (first[i] != second[i]) {
                return false;
            }
        }
        return true;
    }
}

#endif //STRING_STRING_H
