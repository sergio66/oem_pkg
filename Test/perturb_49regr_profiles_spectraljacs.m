chans38 = [        41          54         181         273         317         359         445         449 ...
                  532         758         903         904        1000        1020        1034        1055 ...
                 1075        1103        1249        1282        1291        1447        1475        1557 ...
                 1604        1614        1618        1660        1790        1866        1867        1868 ...
                 1878        1888        2112        2140        2321        2333];
chans41 = [chans38 2325        2339        2353]; chans41 = sort(chans41);

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

[h49,ha49,p49,pa49] = rtpread('/asl/s1/sergio/pin_feb2002_sea_airsnadir_ip.so2.rtp');

rtpwrite('junk.ip.rtp',h49,ha49,p49,pa49);
klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp']; eval(klayerser);
sartaer   = ['!' sarta   ' fin=junk.op.rtp fout=junk.rp.rtp']; eval(sartaer);
[h49x,ha49x,p49x,pa49x] = rtpread('junk.rp.rtp');

mmw = mmwater_rtp(h49x,p49x);
[Y,I] = sort(p49.stemp);
plot(p49.stemp(I),mmw(I),'o-')

p49new = perturb_water(h49,p49,+0.10,-1);
rtpwrite('junk.ip.rtp',h49,ha49,p49new,pa49);
klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp']; eval(klayerser);
sartaer   = ['!' sarta   ' fin=junk.op.rtp fout=junk.rp.rtp']; eval(sartaer);
[h49x,ha49x,p49newx,pa49x] = rtpread('junk.rp.rtp');

error(';kjkj')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_up_oem_pkg_spectraljacs1
set_up_oem_pkg_spectraljacs2
set_up_oem_pkg_spectraljacs_cld

