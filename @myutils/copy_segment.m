function [ns] = copy_segment(s)
%COPY_SEGMENT copying segment.
%   given a segment s, construct a new copy.
    ns = segment();
    
    %copy each property by property
    ns.t         = s.t;
    ns.w         = s.w;
    ns.speed     = s.speed;
    ns.kappa     = s.kappa;
    ns.relaccel  = s.relaccel;
    ns.qtype     = s.qtype;
    ns.eloc      = s.eloc;
    ns.eang      = s.eang;
    ns.Z         = s.Z;
    ns.Zp        = s.Zp;
    ns.Zpp       = s.Zpp;
    ns.Zn        = s.Zn;
    ns.approxv   = s.approxv;
    ns.dom       = s.dom;
    ns.bcside    = s.bcside;
    ns.a         = s.a;
    ns.b         = s.b;
    ns.f         = s.f;
    ns.g         = s.g;
    ns.qpblocha  = s.qpblocha;
end

