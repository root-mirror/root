
namespace Outer
{
    class Inner
    {
    public:
        Inner();
        void fn1(int arg);
    };

    template<typename T, int U>
    class Template
    {
    public:
        Template();
    };

    /**
     * Simple specialisation.
     */
    template<>
    class Template<Inner, 1>
    {
    public:
        Template();
    };

    /**
     * Complex specialisation.
     */
    template<>
    class Template<Template<Inner, 2>, 1>
    {
    public:
        Template();
    };

    /**
     * Simply templated function.
     */
    template<typename T, int U>
    void doit(Template<T,U> arg);

    /**
     * Less-simply templated function.
     */
    template<typename T, int U>
    void doit(Template<Template<T,U>,U> arg);
};