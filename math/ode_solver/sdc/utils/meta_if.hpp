#ifndef META_IF_HPP
#define META_IF_HPP

template<bool conditional, typename Then, typename Else>
struct IF
{
    typedef Then type;       
};
template<typename Then, typename Else>
struct IF<false,Then,Else>
{
    typedef Else type;    
};

template<bool conditional, typename Then, typename Else>
struct IFP : public IF<conditional,Then,Else>::type {};


#endif
