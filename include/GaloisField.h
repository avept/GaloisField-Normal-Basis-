#pragma once

#include <bitset>
#include <string_view>
#include <string>

class GaloisField
{
    // magic constants
    using FieldElement = std::vector<bool>;
    inline static constexpr int64_t fieldDimension = 293;
    inline static constexpr int64_t extendedDimension = fieldDimension * 2;

    static std::unordered_map<std::size_t, std::size_t> positions;
    static std::vector<std::vector<std::size_t>> multiplicativeMatrix;

public:
    GaloisField() = default;
    GaloisField(const std::string_view& element);
    GaloisField(const FieldElement& bits);

    GaloisField operator+(const GaloisField& rhs) const;
    GaloisField operator*(const GaloisField& rhs) const;
    GaloisField& operator =(const GaloisField& rhs) = default;

    GaloisField square() const noexcept;
    bool trace() const noexcept;
    GaloisField power(const GaloisField& degree) const noexcept;
    GaloisField inverse() const noexcept;

    static void createMultiplicativeMatrix() noexcept;

    friend std::ostream& operator <<(std::ostream& out, const GaloisField& rhs);

private:
    // functions
    GaloisField shiftBitsToLow(const int64_t count) const noexcept;

    static int64_t evaluateMatrixElement(const std::size_t i, const std::size_t j) noexcept;
    static int64_t findMultiplicativeElement(const std::size_t position, const int64_t module) noexcept;

    bool evaluateBit(const GaloisField& lhs, const GaloisField& rhs) const noexcept;

    // members
    FieldElement m_data;
};