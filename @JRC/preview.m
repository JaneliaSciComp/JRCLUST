function preview(obj)
    %PREVIEW Display Preview GUI
    obj.startParPool();
    hPreview = jrclust.controllers.curate.PreviewController(obj.hCfg);
    hPreview.preview();
end