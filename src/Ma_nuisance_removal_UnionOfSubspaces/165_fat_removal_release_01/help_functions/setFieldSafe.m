function params = setFieldSafe(params, fieldName, v, vDefault)
    if ~isempty(v)
        params  = setfield(params,fieldName,v);
    else
        params  = setfield(params,fieldName, vDefault);
    end
end