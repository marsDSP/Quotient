#pragma once

#include "PluginProcessor.h"

//==============================================================================
class PluginProcessorEditor final : public AudioProcessorEditor
{
public:
    explicit PluginProcessorEditor (PluginProcessor&);
    ~PluginProcessorEditor() override;

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;

private:
    PluginProcessor& pRef;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PluginProcessorEditor)
};
