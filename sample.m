 function [output]= sample(delta, p)

    index = dsearchn(p.benthic_sp_2, delta);

    md = p.pdfs(index).mynewfield;

    output = random(md,1);

end