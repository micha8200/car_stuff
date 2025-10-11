function LegendItemHitFnc(obj, event)
if strcmp(event.SelectionType, 'normal') && strcmp(event.Region, 'label') && isprop(event, 'Peer') && ~isempty(event.Peer)
    if strcmp(event.Peer.Visible, 'on')
        event.Peer.Visible = 'off';
    else
        event.Peer.Visible = 'on';
    end
end
end