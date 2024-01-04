#include <GaloisField.h>
#include <iostream>
#include <map>
#include <cassert>

std::unordered_map<std::size_t, std::size_t> GaloisField::positions;
std::vector<std::vector<std::size_t>> GaloisField::multiplicativeMatrix;

GaloisField::GaloisField(const std::string_view& element)
{
    if(element.size() != GaloisField::fieldDimension)
    {
        std::cout << "The available length is exceeded!" << std::endl;
        assert(false);
    }

    m_data.resize(GaloisField::fieldDimension, false);
    for(int64_t i = 0; i < element.size(); ++i)
    {
        m_data[i] = (element[i] == '1');
    }
}

GaloisField::GaloisField(const FieldElement& bits) : m_data(bits)
{
    // 
}

GaloisField GaloisField::operator+(const GaloisField& rhs) const
{
    std::vector<bool> result(GaloisField::fieldDimension);
    for(std::size_t i = 0; i < GaloisField::fieldDimension; ++i)
    {
        result[i] = m_data[i] ^ rhs.m_data[i];
    }

    return result;
}

GaloisField GaloisField::operator*(const GaloisField& rhs) const
{
    GaloisField lhsCopy = *this;
    GaloisField rhsCopy = rhs;

    std::vector<bool> result;
    for (std::size_t i = 0; i < GaloisField::fieldDimension; ++i) 
    {
        result.push_back(evaluateBit(lhsCopy, rhsCopy));

        lhsCopy = lhsCopy.shiftBitsToLow(1);
        rhsCopy = rhsCopy.shiftBitsToLow(1);
    }

    return result;
}

std::ostream& operator<<(std::ostream& out, const GaloisField& rhs)
{
    for(std::size_t i = 0; i < 293; ++i)
    {
        out << (rhs.m_data[i] ? '1' : '0');
    }
    return out;
}

bool GaloisField::trace() const noexcept
{
    bool result = false;
    for(const auto& it : m_data)
    {
        result ^= it;
    }

    return result;
}

GaloisField GaloisField::square() const noexcept
{
    std::vector<bool> result(GaloisField::fieldDimension, false);
    result[0] = m_data[GaloisField::fieldDimension - 1];

    for (std::size_t i = GaloisField::fieldDimension - 1; i > 0; --i)
    {
        result[i] = m_data[i - 1];
    }

    return result;
}

GaloisField GaloisField::power(const GaloisField& degree) const noexcept
{
    std::vector<bool> neutralCoeffs(GaloisField::fieldDimension, true);
    GaloisField result(neutralCoeffs);
    GaloisField base = *this;

    if (!degree.m_data.empty() && degree.m_data[0]) 
    {
        result = result * base;
    }

    for (std::size_t i = 1; i < degree.m_data.size(); ++i) 
    {
        result = result.square();
        if (degree.m_data[i]) 
        {
            result = result * base;
        }
    }

    return result;
}

GaloisField GaloisField::inverse() const noexcept
{
    GaloisField beta = *this;

    int k = 1;
    std::string m_binary = "100100100"; // m - 1 = 292

    for (int i = 1; i <= 8; ++i)
    {
        GaloisField original_beta = beta;

        for (int j = 0; j < k; ++j)
        {
            beta = beta.square();
        }

        beta = beta * original_beta;
        k *= 2;

        if (m_binary[i] == '1') 
        {
            GaloisField squared_beta = beta.square();
            beta = squared_beta * (*this);
            ++k;
        }
    }
    
    beta = beta.square();
    return beta;
}

void GaloisField::createMultiplicativeMatrix() noexcept
{
    int64_t numberOfOne = 0;
    multiplicativeMatrix.resize(GaloisField::fieldDimension);

    for (std::size_t i = 0; i < GaloisField::fieldDimension; ++i)
    {
        numberOfOne = 0;
        for (std::size_t j = 0; j < GaloisField::fieldDimension; ++j) 
        {
            if (evaluateMatrixElement(i, j) == 1) 
            {
                multiplicativeMatrix[i].push_back(j);
            }
        }
    }
}

GaloisField GaloisField::shiftBitsToLow(const int64_t count) const noexcept
{
    std::vector<bool> result(m_data.size());
    for (std::size_t i = 0; i < m_data.size(); ++i) 
    {
        std::size_t oldPosition = (i + count) % m_data.size();
        result[i] = m_data[oldPosition];
    }

    return result;
}

int64_t GaloisField::evaluateMatrixElement(const std::size_t i, const std::size_t j) noexcept
{
    int64_t module = GaloisField::extendedDimension + 1;
    int64_t pow_i = findMultiplicativeElement(i, module);
    int64_t pow_j = findMultiplicativeElement(j, module);

    std::vector<int64_t> results = { (pow_i + pow_j) % module,
                                     (pow_i - pow_j + module) % module,
                                     (-pow_i + pow_j + module) % module,
                                     (-pow_i - pow_j + module + module) % module };

    for (const auto& it : results) 
    {
        if (it == 1 || it == (module - 1))
            return 1;
    }

    return 0;
}

int64_t GaloisField::findMultiplicativeElement(const std::size_t position, const int64_t module) noexcept
{
    if (positions.find(position) != positions.end()) 
    {
        return positions[position];
    }

    int64_t result = 1;
    for (std::size_t i = 0; i < position; ++i) 
    {
        result = (result * 2) % module;
    }

    positions[position] = result;
    return result;
}

bool GaloisField::evaluateBit(const GaloisField& lhs, const GaloisField& rhs) const noexcept
{
    bool result = lhs.m_data[0] & rhs.m_data[1]; // hardcode
    for(std::size_t i = 1; i < 293; ++i)
    {
        result ^= lhs.m_data[i] & (rhs.m_data[multiplicativeMatrix[i][0]] ^ rhs.m_data[multiplicativeMatrix[i][1]]);
    }

    return result;
}
