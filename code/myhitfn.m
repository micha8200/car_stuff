function myhitfn(obj, event, varargin)
if strcmpi(event.Peer.Visible, 'on')
    event.Peer.Visible = 'off';
else
    event.Peer.Visible = 'on';
end
end