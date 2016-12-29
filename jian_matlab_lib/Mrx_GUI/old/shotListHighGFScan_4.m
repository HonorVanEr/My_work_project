function [out,pos,fcshift]=shotListHighGFScan_4
% GF LP pos scan on ?
% high fill pressure (2.5/0.9, 3.6 mTorr), GF -50
% FC shift in positive direction for second half

out=168000+[
    423:430 ...
    431:439 ...
    440:448 ...
    449:456 ...
    457:464 ...
    465:472 ...
    473:480 ...
    481:488 ...%
    490:497 ...
    498:505 ...
    ];

pos=[
    ones(1,8)*41.5 ...
    ones(1,9)*39.5 ...
    ones(1,9)*37.5 ...
    ones(1,8)*35.5 ...
    ones(1,8)*33.5 ...
    ones(1,8)*41.5 ...
    ones(1,8)*39.5 ...
    ones(1,8)*37.5 ...
    ones(1,8)*35.5 ...
    ones(1,8)*33.5 ...
    ];

fcshift=[
    ones(1,8)*+3 ...
    ones(1,9)*+3 ...
    ones(1,9)*+3 ...
    ones(1,8)*+3 ...
    ones(1,8)*+3 ...
    ones(1,8)*0 ...
    ones(1,8)*0 ...
    ones(1,8)*0 ...
    ones(1,8)*0 ...
    ones(1,8)*0 ...
    ];
    