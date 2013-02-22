   % Don't need goodind(iQ1), iQ1 will do..., etc.
   % quick ugly fix when not retrieval water
   if length(iQ1) > 10
     wsum = sum(jac(:,goodind(iQ1)),2);
   else
      wsum = 1;
   end
   
   tsum = sum(jac(:,goodind(itemp)),2);

   wsum_max = max(abs(wsum));
   tsum_max = max(abs(tsum));
   for i = 1 : length(qstjacindex)
      qsum_max(i) = max(abs(jac(:,i)));
   end

   % Pick temperature as the standard
   w_mult = tsum_max/wsum_max;
   for i=1:length(qstjacindex)
      q_mult(i) = tsum_max/qsum_max(i);
   end

   bonk=length(qstjacindex)+1 : length(qstjacindex)+numlays;
   % Now apply to Jacobian and modify qrenorm 
   jac(:,bonk) = jac(:,bonk).*w_mult;
   qrenorm(bonk) = qrenorm(bonk).*w_mult;
   for i = 1 : length(qstjacindex)
      jac(:,i)   = jac(:,i).*q_mult(i);
      qrenorm(i) = qrenorm(i).*q_mult(i);
   end
