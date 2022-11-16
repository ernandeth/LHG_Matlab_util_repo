function rotax(axis)

V=get(gca,'View');
V(axis) = get(gco,'Val');
set(gca,'View', V );
