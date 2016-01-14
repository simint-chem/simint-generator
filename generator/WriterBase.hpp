#ifndef WRITERBASE_HPP
#define WRITERBASE_HPP

class WriterBase
{
public:
    virtual bool IsInline(void) const = 0;
    virtual bool IsExternal(void) const = 0;
};

#endif
