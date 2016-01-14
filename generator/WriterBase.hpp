#ifndef WRITERBASE_HPP
#define WRITERBASE_HPP

class WriterBase
{
public:
    virtual bool IsInline(void) const = 0;
    virtual bool IsExternal(void) const = 0;
};

class IsInlineRR
{
public:
    virtual bool IsInline(void) const { return true; }
    virtual bool IsExternal(void) const { return false; }
};

class IsExternalRR
{
public:
    virtual bool IsInline(void) const { return false; }
    virtual bool IsExternal(void) const { return true; }
};

#endif
