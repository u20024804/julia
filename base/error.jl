# pseudo-definitions to show how everything behaves
#
# throw(label, val) = # throw a value to a dynamically enclosing block
#
# function rethrow(val)
#     global current_exception = val
#     throw(current_handler(), current_exception)
# end
#
# rethrow() = rethrow(current_exception)
#
# function throw(val)
#     global catch_backtrace = backtrace()
#     rethrow(val)
# end

## native julia error handling ##

error(s::AbstractString) = throw(ErrorException(s))
error(s...) = throw(ErrorException(string(s...)))

macro unexpected()
    :(error("unexpected branch reached"))
end

rethrow() = ccall(:jl_rethrow, Void, ())::Bottom
rethrow(e) = ccall(:jl_rethrow_other, Void, (Any,), e)::Bottom
backtrace() = ccall(:jl_backtrace_from_here, Array{Ptr{Void},1}, ())
catch_backtrace() = ccall(:jl_get_backtrace, Array{Ptr{Void},1}, ())

## system error handling ##

errno() = ccall(:jl_errno, Cint, ())
errno(e::Integer) = ccall(:jl_set_errno, Void, (Cint,), e)
strerror(e::Integer) = bytestring(ccall(:strerror, Ptr{UInt8}, (Int32,), e))
strerror() = strerror(errno())
systemerror(p, b::Bool) = b ? throw(SystemError(string(p))) : nothing

## assertion functions and macros ##

assert(x) = x ? nothing : throw(AssertionError())
macro assert(ex, msgs...)
    if ex === true
        :(nothing)
    elseif isempty(msgs)
        if isa(ex, Expr)
            :($(esc(ex)) ? nothing : throw(AssertionError($(Expr(:quote,ex)))))
        elseif ex === false
            :(throw(AssertionError(false)))
        else
            :(throw(TypeError(:toplevel, "", Bool, $(esc(ex)))))
        end
    else
        if isa(ex, Expr)
            :($(esc(ex)) ? nothing : throw(AssertionError($(Expr(:quote,ex)), $(esc(msgs[1])))))
        elseif ex === false
            :(throw(AssertionError(false, $(esc(msgs[1])))))
        else
            :(throw(TypeError(:toplevel, "", Bool, $(esc(ex)))))
        end
    end
end
