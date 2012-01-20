#ifndef META_HPP
#define META_HPP

template<bool conditional, typename Then, typename Else>
struct IF
{
    typedef Then type;
};
template<typename Then, typename Else>
struct IF<false, Then, Else>
{
    typedef Else type;
};

template<bool conditional, typename Then, typename Else>
struct IFP : public IF<conditional, Then, Else>::type {};

template<bool condition>
struct StaticAssert;

template<>
struct StaticAssert<true> { enum {WRONG_QUADRATURE_CHOSEN}; };

#define STATIC_ASSERT(CONDITION,MSG) \
    if (StaticAssert<bool(CONDITION)>::MSG) {}

#define STATIC_ASSERT_QUAD(TYPE) \
    STATIC_ASSERT((TYPE::quadrature == gauss_lobatto \
                || TYPE::quadrature == gauss_radau \
                || TYPE::quadrature == clenshaw_curtis), \
                WRONG_TYPE_OF_QUADRATURE)



#endif
